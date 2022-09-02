!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportComputes(iter, convergence, elemcon)
!-----------------------------------------------------------------------------------------------------------!
use arrays_mod,       only: phi_mx, phi_gr, phi_gr_indiv, ww, ww_new, ww_mix,                         &
                            qmx_final, qgr_final, qmx_interp_mm, qmx_interp_mg, qgr_interp, ds_gr_ed, &
                            ds_mx_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, volnp
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

character(40) :: file_name
!-----------------------------------------------------------------------------------------------------------!
adsorbed = .False.

if (export(exportPhiNodal, iter, convergence)) call ExportNodalProfile(phi_mx, phi_gr, numNodes, nodeCoord, volnp)
if (export(exportField, iter, convergence))    call ExportFieldAscii(ww, ww_new, ww_mix)

if (export(exportPhiSmeared, iter, convergence)) then
  ! Planar surfaces
  do mm = 1, 3
    do nn = 1, 2
      if (isDirichletFace(mm,nn)) then
        file_name = ""
        write(file_name,'("o.phi_smear_w",I1,"_",I1)') mm, nn
        call ExportSmearedProfile(planar_cell_of_np(:,mm,nn), cell_vol_planar(:,mm,nn), numNodes, file_name, phi_mx, phi_gr, volnp, lbin, nbin)
      endif
    enddo
  enddo

  ! Spherical nanoparticles
  do mm = 1, numNanopFaces
    file_name = ""
    write(file_name,'("o.phi_smear_np",I1)') mm
    call ExportSmearedProfile(sph_cell_of_np(mm,:), cell_vol_sph(mm,:), numNodes, file_name, phi_mx, phi_gr, volnp, lbin, nbin)
  enddo
endif

if (matrixExist.eq.1) then
  if (export(exportPhiEndMiddle, iter, convergence)) call ExportEndMiddleProfile(numEdwPointsMatrix, qmx_final, qmx_final, "mx", numNodes, nodeCoord)
  if (export(exportPropagators, iter, convergence))  call ExportPropagator(numEdwPointsMatrix, qmx_final, "mx")

  ! Planar surfaces
  do mm = 1, 3
    do nn = 1, 2
      if (isDirichletFace(mm,nn)) then
        if (export(exportChainsPerArea, iter, convergence)) then
          file_name = ""
          write(file_name,'("o.chains_area_w",I1,"_",I1)') mm, nn
          call ExportChainsArea(nodeBelongsToDirichletFace, elemcon, planar_cell_of_np(:,mm,nn), "mx", rg2OfMatrixMonomer, lengthMatrix, numEdwPointsMatrix, ds_mx_ed, qmx_final, phi_mx, ww)
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
      file_name = ""
      write(file_name,'("o.chains_area_w",I1,"_",I1)') mm, nn
      call ExportChainsArea(nodeBelongsToDirichletFace, elemcon, sph_cell_of_np(mm,:), "mx", rg2OfMatrixMonomer, lengthMatrix, numEdwPointsMatrix, ds_mx_ed, qmx_final, phi_mx, ww)
    endif
    if (export(exportAdsorbedFree, iter, convergence)) then
      do kk = 1, numNodes
        if (dist_from_np(mm,kk)<adsorptionDistance) adsorbed(kk) = .true.
      enddo
    endif
  enddo

  if (export(exportAdsorbedFree, iter, convergence)) call ExportAdsorbed(nodeBelongsToDirichletFace, elemcon, adsorbed)

#ifdef DEBUG_OUTPUTS
  call ExportPropagator(numConvolPointsMatrix, qmx_interp_mm, "mm")
  call ExportPropagator(numConvolPointsGrafted, qmx_interp_mg, "mg")
  write(6,'(2X,A40)')adjl("****************************************",40)
#endif
endif

if (graftedExist.eq.1) then
  if (export(exportPhiEndMiddle, iter, convergence)) call ExportEndMiddleProfile(numConvolPointsGrafted, qgr_interp, qmx_interp_mg, "gr", numNodes, nodeCoord)
  if (export(exportPropagators, iter, convergence))  call ExportPropagator(numEdwPointsGrafted, qgr_final, "gr")

  if (export(exportPhiIndividual, iter, convergence)) then
    call ComputeIndivProfile(numNodes, elemcon, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, ww, targetNumGraftedChains, graftPointId, graftPointValue, phi_gr_indiv)

    if (exportAllGraftedChains.eq.1) then
      call ExportIndivProfile(targetNumGraftedChains, numNodes, nodeCoord, phi_gr_indiv)
    else
      call ExportIndivProfile(numGraftedChainsToExport, numNodes, nodeCoord, phi_gr_indiv)
    endif

    call ExportVtuIndiv(numGraftedChainsToExport, phi_gr_indiv)
  endif

  ! Planar surfaces
  do mm = 1, 3
    do nn = 1, 2
      if (isDirichletFace(mm,nn)) then
        if (export(exportBrushThickness, iter, convergence)) then
          file_name = ""
          write(file_name,'("o.brush_w",I1,"_",I1)') mm, nn
          call ExportBrush(targetNumGraftedChains, numNodes, phi_gr, phi_gr_indiv, volnp, file_name, dist_from_face(:,mm,nn))
          file_name = ""
          write(file_name,'("o.brush99_w",I1,"_",I1)') mm, nn
          call ExportBrush99(planar_cell_of_np(:,mm,nn), targetNumGraftedChains, numNodes, file_name, phi_gr, phi_gr_indiv, volnp, lbin, nbin)
        endif
        if (export(exportChainsPerArea, iter, convergence)) then
          file_name = ""
          write(file_name,'("o.chains_area_w",I1,"_",I1)') mm, nn
          call ExportChainsArea(nodeBelongsToDirichletFace, elemcon, planar_cell_of_np(:,mm,nn), "gr", rg2OfGraftedMonomer, lengthGrafted, numEdwPointsGrafted, ds_gr_ed, qgr_final, phi_gr, ww)
        endif
      endif
    enddo
  enddo

  ! Spherical nanoparticles
  do mm = 1, numNanopFaces
    if (export(exportBrushThickness, iter, convergence)) then
      file_name = ""
      write(file_name,'("o.brush_np",I1)') mm
      call ExportBrush(targetNumGraftedChains, numNodes, phi_gr, phi_gr_indiv, volnp, file_name, dist_from_np(mm,:))
      file_name = ""
      write(file_name,'("o.brush99_np",I1)') mm
      call ExportBrush99(sph_cell_of_np(mm,:), targetNumGraftedChains, numNodes, file_name, phi_gr, phi_gr_indiv, volnp, lbin, nbin)
    endif
    if (export(exportChainsPerArea, iter, convergence)) then
      file_name = ""
      write(file_name,'("o.chains_area_w",I1,"_",I1)') mm, nn
      call ExportChainsArea(nodeBelongsToDirichletFace, elemcon, sph_cell_of_np(mm,:), "gr", rg2OfGraftedMonomer, lengthGrafted, numEdwPointsGrafted, ds_gr_ed, qgr_final, phi_gr, ww)
    endif
  enddo
#ifdef DEBUG_OUTPUTS
  call ExportPropagator(numConvolPointsGrafted, qgr_interp, "gg")
#endif
endif

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine ExportComputes
