!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportComputes(iter, convergence, elemcon)
!-----------------------------------------------------------------------------------------------------------!
use arrays_mod,       only: phiMatrix, phiGrafted, phiGraftedIndiv,                     &
                            wwField, wwFieldNew, wwFieldMixed,                          &
                            qqMatrixFinal, qqMatrixInterp, qqMatrixInterpGrafted,       &
                            qqGraftedFinal, qqGraftedInterp, dsEdwGrafted,              &
                            dsEdwMatrix, xsEdwGrafted, xsConvGrafted, coeffConvGrafted, &
                            nodeVolume
use hist_mod,         only: numBins, binLength, planarCellId, sphericalCellId, distanceFromFace, &
                            distanceFromNanop, planarCellVolume, sphericalCellVolume
use delta_mod,        only: graftPointValue, graftPointId, targetNumGraftedChains
use geometry_mod,     only: numNodes, nodeCoord, isDirichletFace, nodeBelongsToDirichletFace
use write_helper_mod, only: adjl, export
use parser_vars_mod,  only: numNanopFaces, matrixExist, graftedExist, adsorptionDistance,         &
                            rg2OfMatrixMonomer, rg2OfGraftedMonomer, lengthMatrix, lengthGrafted, &
                            numEdwPointsMatrix, numConvolPointsMatrix,                            &
                            numEdwPointsGrafted, numConvolPointsGrafted,                          &
                            exportPhiNodal,                                                       &
                            exportPhiSmeared,                                                     &
                            exportPhiEndMiddle,                                                   &
                            exportField,                                                          &
                            exportPropagators,                                                    &
                            exportPhiIndividual,                                                  &
                            exportBrushThickness,                                                 &
                            exportChainsPerArea,                                                  &
                            exportAdsorbedFree,                                                   &
                            exportAllGraftedChains,                                               &
                            numGraftedChainsToExport
use fhash_module__ints_double
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
logical, intent(in)          :: convergence
logical, dimension(numNodes) :: adsorbed

integer, intent(in) :: iter
integer             :: node, axis, face

type(fhash_type__ints_double), intent(inout) :: elemcon

character(40) :: fileName
!-----------------------------------------------------------------------------------------------------------!
adsorbed = .False.

if (export(exportPhiNodal, iter, convergence)) CALL ExportNodalProfile(phiMatrix, phiGrafted, numNodes, nodeCoord, nodeVolume)
if (export(exportField, iter, convergence))    CALL ExportFieldAscii(wwField, wwFieldNew, wwFieldMixed)

if (export(exportPhiSmeared, iter, convergence)) then
  ! Planar surfaces
  do axis = 1, 3
    do face = 1, 2
      if (isDirichletFace(axis,face)) then
        fileName = ""
        write(fileName,'("o.phi_smear_w",I1,"_",I1)') axis, face
        CALL ExportSmearedProfile(planarCellId(:,axis,face), planarCellVolume(:,axis,face), numNodes, fileName, phiMatrix, phiGrafted, nodeVolume, binLength, numBins)
      endif
    enddo
  enddo

  ! Spherical nanoparticles
  do face = 1, numNanopFaces
    fileName = ""
    write(fileName,'("o.phi_smear_np",I1)') face
    CALL ExportSmearedProfile(sphericalCellId(face,:), sphericalCellVolume(face,:), numNodes, fileName, phiMatrix, phiGrafted, nodeVolume, binLength, numBins)
  enddo
endif

if (matrixExist.eq.1) then
  if (export(exportPhiEndMiddle, iter, convergence)) CALL ExportEndMiddleProfile(numEdwPointsMatrix, qqMatrixFinal, qqMatrixFinal, "mx", numNodes, nodeCoord)
  if (export(exportPropagators, iter, convergence))  CALL ExportPropagator(numEdwPointsMatrix, qqMatrixFinal, "mx")

  ! Planar surfaces
  do axis = 1, 3
    do face = 1, 2
      if (isDirichletFace(axis,face)) then
        if (export(exportChainsPerArea, iter, convergence)) then
          fileName = ""
          write(fileName,'("o.chains_area_w",I1,"_",I1)') axis, face
          CALL ExportChainsArea(nodeBelongsToDirichletFace, elemcon, planarCellId(:,axis,face), "mx", rg2OfMatrixMonomer, &
                                lengthMatrix, numEdwPointsMatrix, dsEdwMatrix, qqMatrixFinal, phiMatrix, wwField)
        endif
        if (export(exportAdsorbedFree, iter, convergence)) then
          do node = 1, numNodes
            if (distanceFromFace(node,axis,face)<adsorptionDistance) adsorbed(node) = .True.
          enddo
        endif
      endif
    enddo
  enddo

  ! Spherical nanoparticles
  do face = 1, numNanopFaces
    if (export(exportChainsPerArea, iter, convergence)) then
      fileName = ""
      write(fileName,'("o.chains_area_w",I1,"_",I1)') face
      CALL ExportChainsArea(nodeBelongsToDirichletFace, elemcon, sphericalCellId(face,:), "mx", rg2OfMatrixMonomer, lengthMatrix, &
                            numEdwPointsMatrix, dsEdwMatrix, qqMatrixFinal, phiMatrix, wwField)
    endif
    if (export(exportAdsorbedFree, iter, convergence)) then
      do node = 1, numNodes
        if (distanceFromNanop(face,node) < adsorptionDistance) adsorbed(node) = .True.
      enddo
    endif
  enddo

  if (export(exportAdsorbedFree, iter, convergence)) CALL ExportAdsorbed(nodeBelongsToDirichletFace, elemcon, adsorbed)

#ifdef DEBUG_OUTPUTS
  CALL ExportPropagator(numConvolPointsMatrix, qqMatrixInterp, "mm")
  CALL ExportPropagator(numConvolPointsGrafted, qqMatrixInterpGrafted, "mg")
  write(6,'(2X,A40)')adjl("****************************************",40)
#endif
endif

if (graftedExist.eq.1) then
  if (export(exportPhiEndMiddle, iter, convergence)) CALL ExportEndMiddleProfile(numConvolPointsGrafted, qqGraftedInterp, qqMatrixInterpGrafted, "gr", numNodes, nodeCoord)
  if (export(exportPropagators, iter, convergence))  CALL ExportPropagator(numEdwPointsGrafted, qqGraftedFinal, "gr")

  if (export(exportPhiIndividual, iter, convergence)) then
    CALL ComputeIndivProfile(numNodes, elemcon, qqMatrixInterpGrafted, dsEdwGrafted, xsEdwGrafted, xsConvGrafted, coeffConvGrafted, wwField, targetNumGraftedChains, graftPointId, graftPointValue, phiGraftedIndiv)

    if (exportAllGraftedChains.eq.1) then
      CALL ExportIndivProfile(targetNumGraftedChains, numNodes, nodeCoord, phiGraftedIndiv)
    elseif ((exportAllGraftedChains.eq.0).or.(exportAllGraftedChains.eq.2)) then
      CALL ExportIndivProfile(numGraftedChainsToExport, numNodes, nodeCoord, phiGraftedIndiv)
    endif

    CALL ExportVtuIndiv(numGraftedChainsToExport, phiGraftedIndiv)
  endif

  ! Planar surfaces
  do axis = 1, 3
    do face = 1, 2
      if (isDirichletFace(axis,face)) then
        if (export(exportBrushThickness, iter, convergence)) then
          fileName = ""
          write(fileName,'("o.brush_w",I1,"_",I1)') axis, face
          CALL ExportBrush(targetNumGraftedChains, numNodes, phiGrafted, phiGraftedIndiv, nodeVolume, fileName, distanceFromFace(:,axis,face))
          fileName = ""
          write(fileName,'("o.brush99_w",I1,"_",I1)') axis, face
          CALL ExportBrush99(planarCellId(:,axis,face), targetNumGraftedChains, numNodes, fileName, phiGrafted, phiGraftedIndiv, nodeVolume, binLength, numBins)
        endif
        if (export(exportChainsPerArea, iter, convergence)) then
          fileName = ""
          write(fileName,'("o.chains_area_w",I1,"_",I1)') axis, face
          CALL ExportChainsArea(nodeBelongsToDirichletFace, elemcon, planarCellId(:,axis,face), "gr", rg2OfGraftedMonomer, &
                                lengthGrafted, numEdwPointsGrafted, dsEdwGrafted, qqGraftedFinal, phiGrafted, wwField)
        endif
      endif
    enddo
  enddo

  ! Spherical nanoparticles
  do face = 1, numNanopFaces
    if (export(exportBrushThickness, iter, convergence)) then
      fileName = ""
      write(fileName,'("o.brush_np",I1)') face
      CALL ExportBrush(targetNumGraftedChains, numNodes, phiGrafted, phiGraftedIndiv, nodeVolume, fileName, distanceFromNanop(face,:))
      fileName = ""
      write(fileName,'("o.brush99_np",I1)') face
      CALL ExportBrush99(sphericalCellId(face,:), targetNumGraftedChains, numNodes, fileName, phiGrafted, phiGraftedIndiv, nodeVolume, binLength, numBins)
    endif
    if (export(exportChainsPerArea, iter, convergence)) then
      fileName = ""
      write(fileName,'("o.chains_area_w",I1,"_",I1)') face
      CALL ExportChainsArea(nodeBelongsToDirichletFace, elemcon, sphericalCellId(face,:), "gr", rg2OfGraftedMonomer, lengthGrafted, &
                            numEdwPointsGrafted, dsEdwGrafted, qqGraftedFinal, phiGrafted, wwField)
    endif
  enddo
#ifdef DEBUG_OUTPUTS
  CALL ExportPropagator(numConvolPointsGrafted, qqGraftedInterp, "gg")
#endif
endif

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine ExportComputes
