!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ComputeIndivProfile(numNodes, elemcon, qqMatrixInterpGrafted, dsEdwGrafted, xsEdwGrafted, xsConvGrafted, coeffConvGrafted, ww, &
                             targetNumGraftedChains, graftingPointId, graftingPointInitialValue, phiGraftedIndiv)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: nodeBelongsToDirichletFace
use write_helper_mod, only: adjl
use parser_vars_mod,  only: numConvolPointsGrafted, numEdwPointsGrafted, lengthGrafted, &
                            rg2OfGraftedMonomer, exportAllGraftedChains,                &
                            numGraftedChainsToExport, graftingPointIndexToExport
use geometry_mod,     only: nodeBelongsToDirichletFace
use write_helper_mod, only: adjl
use error_handling_mod
use fhash_module__ints_double
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                                    :: numNodes, targetNumGraftedChains
integer, intent(in), dimension(targetNumGraftedChains) :: graftingPointId
integer                                                :: graftingPointIndex, graftedChain, node

type(fhash_type__ints_double), intent(inout) :: elemcon

real(8), intent(in), dimension(targetNumGraftedChains)            :: graftingPointInitialValue
real(8), intent(in), dimension(numNodes)                          :: ww
real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qqMatrixInterpGrafted
real(8), intent(in), dimension(numEdwPointsGrafted+1)             :: dsEdwGrafted, xsEdwGrafted
real(8), intent(in), dimension(numConvolPointsGrafted+1)          :: xsConvGrafted, coeffConvGrafted
real(8), intent(out), dimension(numNodes,targetNumGraftedChains)  :: phiGraftedIndiv
real(8), dimension(2,numNodes)                                    :: qqGrafted
real(8), dimension(numEdwPointsGrafted+1,numNodes)                :: qqGraftedFinal
real(8), dimension(numConvolPointsGrafted+1,numNodes)             :: qqGraftedInterp
!------------------------------------------------------------------------------------------------------!
write(6,'(2X,A43)')adjl("Computing indiv profiles of grafted chains.",43)

CALL FemMatrixAssemble(rg2OfGraftedMonomer, ww)

if (exportAllGraftedChains.eq.1) then ! Export profile of all grafted chains separately
  do graftedChain = 1, targetNumGraftedChains
    qqGrafted      = 0.0d0
    qqGraftedFinal = 0.0d0

    qqGrafted(1,graftingPointId(graftedChain))      = graftingPointInitialValue(graftedChain)
    qqGraftedFinal(1,graftingPointId(graftedChain)) = graftingPointInitialValue(graftedChain)

    write(6, '(2X,A21,1X,I3,1X,A8,1X,I7,1X,A2,1X)', advance='no') "Grafting point index:", graftedChain, "with id:", graftingPointId(graftedChain), "->"
    CALL SolverEdwards(dsEdwGrafted, numEdwPointsGrafted, qqGrafted, qqGraftedFinal, nodeBelongsToDirichletFace, elemcon)

    do node = 1, numNodes
      CALL interp_linear(1, numEdwPointsGrafted+1, xsEdwGrafted, qqGraftedFinal(:,node), numConvolPointsGrafted+1, xsConvGrafted, qqGraftedInterp(:,node))
    enddo

    CALL ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeffConvGrafted, qqGraftedInterp, qqMatrixInterpGrafted, phiGraftedIndiv(:,graftedChain))
  enddo
elseif (exportAllGraftedChains.eq.0) then ! Export profile of specific grafted chains separately
  do graftedChain = 1, numGraftedChainsToExport
    qqGrafted      = 0.0d0
    qqGraftedFinal = 0.0d0

    graftingPointIndex = graftingPointIndexToExport(graftedChain)

    qqGrafted(1,graftingPointId(graftingPointIndex))      = graftingPointInitialValue(graftingPointIndex)
    qqGraftedFinal(1,graftingPointId(graftingPointIndex)) = graftingPointInitialValue(graftingPointIndex)

    write(6, '(2X,A21,1X,I3,1X,A8,1X,I7,1X,A2,1X)', advance='no') "Grafting point index:", graftingPointIndex, "with id:", graftingPointId(graftingPointIndex), "->"
    CALL SolverEdwards(dsEdwGrafted, numEdwPointsGrafted, qqGrafted, qqGraftedFinal, nodeBelongsToDirichletFace, elemcon)

    do node = 1, numNodes
      CALL interp_linear(1, numEdwPointsGrafted+1, xsEdwGrafted, qqGraftedFinal(:,node), numConvolPointsGrafted+1, xsConvGrafted, qqGraftedInterp(:,node))
    enddo

    CALL ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeffConvGrafted, qqGraftedInterp, qqMatrixInterpGrafted, phiGraftedIndiv(:,graftingPointIndex))
  enddo
elseif (exportAllGraftedChains.eq.2) then ! Export profile of specific grafted chains as a sum
  qqGrafted      = 0.0d0
  qqGraftedFinal = 0.0d0

  do graftedChain = 1, numGraftedChainsToExport
    graftingPointIndex = graftingPointIndexToExport(graftedChain)

    qqGrafted(1,graftingPointId(graftingPointIndex))      = graftingPointInitialValue(graftingPointIndex)
    qqGraftedFinal(1,graftingPointId(graftingPointIndex)) = graftingPointInitialValue(graftingPointIndex)
  enddo

  write(6, '(2X,A40)') "Sum of specific grafted chains:"
  CALL SolverEdwards(dsEdwGrafted, numEdwPointsGrafted, qqGrafted, qqGraftedFinal, nodeBelongsToDirichletFace, elemcon)

  do node = 1, numNodes
    CALL interp_linear(1, numEdwPointsGrafted+1, xsEdwGrafted, qqGraftedFinal(:,node), numConvolPointsGrafted+1, xsConvGrafted, qqGraftedInterp(:,node))
  enddo

  CALL ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeffConvGrafted, qqGraftedInterp, qqMatrixInterpGrafted, phiGraftedIndiv(:,1))
endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine ComputeIndivProfile
