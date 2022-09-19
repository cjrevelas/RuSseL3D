!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportComputes(iter, convergence, elemcon)
!-----------------------------------------------------------------------------------------------------------!
use arrays_mod,       only: phiMatrix, phiGrafted, phiGraftedIndiv,               &
                            wwField, wwFieldNew, wwFieldMixed,                    &
                            qqMatrixFinal, qqMatrixInterp, qqMatrixInterpGrafted, &
                            qqGraftedFinal, qqGraftedInterp, ds_gr_ed,            &
                            ds_mx_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv,        &
                            nodeVolume
use hist_mod,         only: nbin, lbin, planar_cell_of_np, sph_cell_of_np, dist_from_face, &
                            dist_from_np, cell_vol_planar, cell_vol_sph
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
integer             :: kk, mm, nn

type(fhash_type__ints_double), intent(inout) :: elemcon

character(40) :: fileName
!-----------------------------------------------------------------------------------------------------------!
adsorbed = .False.

if (export(exportPhiNodal, iter, convergence)) call ExportNodalProfile(phiMatrix, phiGrafted, numNodes, nodeCoord, nodeVolume)
if (export(exportField, iter, convergence))    call ExportFieldAscii(wwField, wwFieldNew, wwFieldMixed)

if (export(exportPhiSmeared, iter, convergence)) then
  ! Planar surfaces
  do mm = 1, 3
    do nn = 1, 2
      if (isDirichletFace(mm,nn)) then
        fileName = ""
        write(fileName,'("o.phi_smear_w",I1,"_",I1)') mm, nn
        call ExportSmearedProfile(planar_cell_of_np(:,mm,nn), cell_vol_planar(:,mm,nn), numNodes, fileName, phiMatrix, phiGrafted, nodeVolume, lbin, nbin)
      endif
    enddo
  enddo

  ! Spherical nanoparticles
  do mm = 1, numNanopFaces
    fileName = ""
    write(fileName,'("o.phi_smear_np",I1)') mm
    call ExportSmearedProfile(sph_cell_of_np(mm,:), cell_vol_sph(mm,:), numNodes, fileName, phiMatrix, phiGrafted, nodeVolume, lbin, nbin)
  enddo
endif

if (matrixExist.eq.1) then
  if (export(exportPhiEndMiddle, iter, convergence)) call ExportEndMiddleProfile(numEdwPointsMatrix, qqMatrixFinal, qqMatrixFinal, "mx", numNodes, nodeCoord)
  if (export(exportPropagators, iter, convergence))  call ExportPropagator(numEdwPointsMatrix, qqMatrixFinal, "mx")

  ! Planar surfaces
  do mm = 1, 3
    do nn = 1, 2
      if (isDirichletFace(mm,nn)) then
        if (export(exportChainsPerArea, iter, convergence)) then
          fileName = ""
          write(fileName,'("o.chains_area_w",I1,"_",I1)') mm, nn
          call ExportChainsArea(nodeBelongsToDirichletFace, elemcon, planar_cell_of_np(:,mm,nn), "mx", rg2OfMatrixMonomer, &
          lengthMatrix, numEdwPointsMatrix, ds_mx_ed, qqMatrixFinal, phiMatrix, wwField)
        endif
        if (export(exportAdsorbedFree, iter, convergence)) then
          do kk = 1, numNodes
            if (dist_from_face(kk,mm,nn)<adsorptionDistance) adsorbed(kk) = .true.
          enddo
        endif
      endif
    enddo
  enddo

  ! Spherical nanoparticles
  do mm = 1, numNanopFaces
    if (export(exportChainsPerArea, iter, convergence)) then
      fileName = ""
      write(fileName,'("o.chains_area_w",I1,"_",I1)') mm, nn
      call ExportChainsArea(nodeBelongsToDirichletFace, elemcon, sph_cell_of_np(mm,:), "mx", rg2OfMatrixMonomer, lengthMatrix, &
                            numEdwPointsMatrix, ds_mx_ed, qqMatrixFinal, phiMatrix, wwField)
    endif
    if (export(exportAdsorbedFree, iter, convergence)) then
      do kk = 1, numNodes
        if (dist_from_np(mm,kk)<adsorptionDistance) adsorbed(kk) = .true.
      enddo
    endif
  enddo

  if (export(exportAdsorbedFree, iter, convergence)) call ExportAdsorbed(nodeBelongsToDirichletFace, elemcon, adsorbed)

#ifdef DEBUG_OUTPUTS
  call ExportPropagator(numConvolPointsMatrix, qqMatrixInterp, "mm")
  call ExportPropagator(numConvolPointsGrafted, qqMatrixInterpGrafted, "mg")
  write(6,'(2X,A40)')adjl("****************************************",40)
#endif
endif

if (graftedExist.eq.1) then
  if (export(exportPhiEndMiddle, iter, convergence)) call ExportEndMiddleProfile(numConvolPointsGrafted, qqGraftedInterp, qqMatrixInterpGrafted, "gr", numNodes, nodeCoord)
  if (export(exportPropagators, iter, convergence))  call ExportPropagator(numEdwPointsGrafted, qqGraftedFinal, "gr")

  if (export(exportPhiIndividual, iter, convergence)) then
    call ComputeIndivProfile(numNodes, elemcon, qqMatrixInterpGrafted, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, wwField, targetNumGraftedChains, graftPointId, graftPointValue, phiGraftedIndiv)

    if (exportAllGraftedChains.eq.1) then
      call ExportIndivProfile(targetNumGraftedChains, numNodes, nodeCoord, phiGraftedIndiv)
    else
      call ExportIndivProfile(numGraftedChainsToExport, numNodes, nodeCoord, phiGraftedIndiv)
    endif

    call ExportVtuIndiv(numGraftedChainsToExport, phiGraftedIndiv)
  endif

  ! Planar surfaces
  do mm = 1, 3
    do nn = 1, 2
      if (isDirichletFace(mm,nn)) then
        if (export(exportBrushThickness, iter, convergence)) then
          fileName = ""
          write(fileName,'("o.brush_w",I1,"_",I1)') mm, nn
          call ExportBrush(targetNumGraftedChains, numNodes, phiGrafted, phiGraftedIndiv, nodeVolume, fileName, dist_from_face(:,mm,nn))
          fileName = ""
          write(fileName,'("o.brush99_w",I1,"_",I1)') mm, nn
          call ExportBrush99(planar_cell_of_np(:,mm,nn), targetNumGraftedChains, numNodes, fileName, phiGrafted, phiGraftedIndiv, nodeVolume, lbin, nbin)
        endif
        if (export(exportChainsPerArea, iter, convergence)) then
          fileName = ""
          write(fileName,'("o.chains_area_w",I1,"_",I1)') mm, nn
          call ExportChainsArea(nodeBelongsToDirichletFace, elemcon, planar_cell_of_np(:,mm,nn), "gr", rg2OfGraftedMonomer, &
          lengthGrafted, numEdwPointsGrafted, ds_gr_ed, qqGraftedFinal, phiGrafted, wwField)
        endif
      endif
    enddo
  enddo

  ! Spherical nanoparticles
  do mm = 1, numNanopFaces
    if (export(exportBrushThickness, iter, convergence)) then
      fileName = ""
      write(fileName,'("o.brush_np",I1)') mm
      call ExportBrush(targetNumGraftedChains, numNodes, phiGrafted, phiGraftedIndiv, nodeVolume, fileName, dist_from_np(mm,:))
      fileName = ""
      write(fileName,'("o.brush99_np",I1)') mm
      call ExportBrush99(sph_cell_of_np(mm,:), targetNumGraftedChains, numNodes, fileName, phiGrafted, phiGraftedIndiv, nodeVolume, lbin, nbin)
    endif
    if (export(exportChainsPerArea, iter, convergence)) then
      fileName = ""
      write(fileName,'("o.chains_area_w",I1,"_",I1)') mm, nn
      call ExportChainsArea(nodeBelongsToDirichletFace, elemcon, sph_cell_of_np(mm,:), "gr", rg2OfGraftedMonomer, lengthGrafted, &
                            numEdwPointsGrafted, ds_gr_ed, qqGraftedFinal, phiGrafted, wwField)
    endif
  enddo
#ifdef DEBUG_OUTPUTS
  call ExportPropagator(numConvolPointsGrafted, qqGraftedInterp, "gg")
#endif
endif

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine ExportComputes
