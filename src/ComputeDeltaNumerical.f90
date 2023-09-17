!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ComputeDeltaNumerical(elemcon, qqMatrixInterpGrafted, ww, deltaNumerical)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: numNodes, nodeBelongsToDirichletFace
use parser_vars_mod,  only: numConvolPointsGrafted, numEdwPointsGrafted, lengthGrafted, &
                            rg2OfGraftedMonomer, molarBulkDensity
use arrays_mod,       only: dsEdwGrafted, xsEdwGrafted, xsConvGrafted, coeffConvGrafted, nodeVolume
use constants_mod,    only: n_avog, m3_to_A3
use delta_mod,        only: graftPointId, targetNumGraftedChains
use write_helper_mod, only: adjl
use error_handling_mod
use fhash_module__ints_double
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: graftedChain, node

type(fhash_type__ints_double), intent(inout) :: elemcon

real(8), intent(in), dimension(numNodes)                          :: ww
real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qqMatrixInterpGrafted
real(8), intent(out), dimension(targetNumGraftedChains)           :: deltaNumerical
real(8), dimension(targetNumGraftedChains)                        :: deltaAnalytic
real(8), dimension(numNodes)                                      :: phiGrafted
real(8), dimension(2,numNodes)                                    :: qqGrafted
real(8), dimension(numEdwPointsGrafted+1,numNodes)                :: qqGraftedFinal
real(8), dimension(numConvolPointsGrafted+1,numNodes)             :: qqGraftedInterp
real(8)                                                           :: numChainsGrafted   = 0.0d0
real(8)                                                           :: initValueTentative = 0.0d0
!------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("****************************************",40)
write(6,'(2X,A40)')adjl("Updating delta of grafted chains:",40)

deltaAnalytic = 0.0d0

! Analytic delta calculation
do graftedChain = 1, targetNumGraftedChains
  deltaAnalytic(graftedChain) = 1.0d0 / nodeVolume(graftPointId(graftedChain)) * m3_to_A3
enddo

! Numerical delta calculation (i.e., through solution of the Edwards equation)
CALL FemMatrixAssemble(rg2OfGraftedMonomer, ww)

do graftedChain = 1, targetNumGraftedChains
  qqGrafted      = 0.0d0
  qqGraftedFinal = 0.0d0

  initValueTentative = deltaAnalytic(graftedChain) * lengthGrafted * 1.0d0 / &
                       (qqMatrixInterpGrafted(numConvolPointsGrafted+1,graftPointId(graftedChain)) * (molarBulkDensity * n_avog))

  qqGrafted(1,graftPointId(graftedChain))      = initValueTentative
  qqGraftedFinal(1,graftPointId(graftedChain)) = initValueTentative

  write(6, '(2X,A19,I7,A3)', advance='no') "Grafting point id: ", graftPointId(graftedChain), " ->"
  CALL SolverEdwards(dsEdwGrafted, numEdwPointsGrafted, qqGrafted, qqGraftedFinal, nodeBelongsToDirichletFace, elemcon)

  do node = 1, numNodes
    CALL interp_linear(1, numEdwPointsGrafted+1, xsEdwGrafted, qqGraftedFinal(:,node), numConvolPointsGrafted+1, xsConvGrafted, qqGraftedInterp(:,node))
  enddo

  CALL ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeffConvGrafted, qqGraftedInterp, qqMatrixInterpGrafted, phiGrafted)

  CALL ComputeNumberOfChains(lengthGrafted, phiGrafted, numChainsGrafted)

  deltaNumerical(graftedChain) = deltaAnalytic(graftedChain) / numChainsGrafted
enddo
write(6,'(2X,A40)')adjl("****************************************",40)

return
!------------------------------------------------------------------------------------------------------!
end subroutine ComputeDeltaNumerical
