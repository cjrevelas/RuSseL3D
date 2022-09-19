!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ComputeDeltaNumerical(elemcon, qmx_interp_mg, ww, deltaNumerical)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: numNodes, nodeBelongsToDirichletFace
use parser_vars_mod,  only: numConvolPointsGrafted, numEdwPointsGrafted, lengthGrafted, &
                            rg2OfGraftedMonomer, molarBulkDensity
use arrays_mod,       only: ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, nodeVolume
use constants_mod,    only: n_avog, m3_to_A3
use delta_mod,        only: graftPointId, targetNumGraftedChains
use write_helper_mod, only: adjl
use error_handing_mod
use fhash_module__ints_double
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: ii, jj

type(fhash_type__ints_double), intent(inout) :: elemcon

real(8), intent(in), dimension(numNodes)                          :: ww
real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qmx_interp_mg
real(8), intent(out), dimension(targetNumGraftedChains)           :: deltaNumerical
real(8), dimension(targetNumGraftedChains)                        :: deltaAnalytic
real(8), dimension(numNodes)                                      :: phiGrafted
real(8), dimension(2,numNodes)                                    :: qgr
real(8), dimension(numEdwPointsGrafted+1,numNodes)                :: qgr_final
real(8), dimension(numConvolPointsGrafted+1,numNodes)             :: qgr_interp
real(8)                                                           :: numChainsGrafted = 0.0d0
real(8)                                                           :: initValueTentative = 0.0d0
!------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("****************************************",40)
write(6,'(2X,A40)')adjl("Updating delta of grafted chains:",40)

deltaAnalytic = 0.0d0

! Analytic delta calculation
do ii = 1, targetNumGraftedChains
  deltaAnalytic(ii) = 1.0d0 / nodeVolume(graftPointId(ii)) * m3_to_A3
enddo

! Numerical delta calculation (i.e., through solution of the Edwards equation)
call FemMatrixAssemble(rg2OfGraftedMonomer, ww)

do ii = 1, targetNumGraftedChains
  qgr       = 0.0d0
  qgr_final = 0.0d0

  initValueTentative = deltaAnalytic(ii) * lengthGrafted * 1.0d0 / &
                       (qmx_interp_mg(numConvolPointsGrafted+1,graftPointId(ii)) * (molarBulkDensity * n_avog))

  qgr(1,graftPointId(ii))       = initValueTentative
  qgr_final(1,graftPointId(ii)) = initValueTentative

  write(6, '(2X,A19,I7,A3)', advance='no') "Grafting point id: ", graftPointId(ii), " ->"
  call SolverEdwards(ds_gr_ed, numEdwPointsGrafted, qgr, qgr_final, nodeBelongsToDirichletFace, elemcon)

  do jj = 1, numNodes
    call interp_linear(1, numEdwPointsGrafted+1, xs_gr_ed, qgr_final(:,jj), numConvolPointsGrafted+1, xs_gr_conv, qgr_interp(:,jj))
  enddo

  call ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeff_gr_conv, qgr_interp, qmx_interp_mg, phiGrafted)

  call ComputeNumberOfChains(lengthGrafted, phiGrafted, numChainsGrafted)

  deltaNumerical(ii) = deltaAnalytic(ii) / numChainsGrafted
enddo
write(6,'(2X,A40)')adjl("****************************************",40)

return
!------------------------------------------------------------------------------------------------------!
end subroutine ComputeDeltaNumerical
