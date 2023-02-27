!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ComputeIndivProfile(numNodes, elemcon, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, ww, &
                             targetNumGraftedChains, gpid, gp_init_value, phi_gr_indiv)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: nodeBelongsToDirichletFace
use write_helper_mod, only: adjl
use parser_vars_mod,  only: numConvolPointsGrafted, numEdwPointsGrafted, lengthGrafted, &
                            rg2OfGraftedMonomer, exportAllGraftedChains,                &
                            numGraftedChainsToExport, gpIndexToExport
use geometry_mod,     only: nodeBelongsToDirichletFace
use write_helper_mod, only: adjl
use error_handing_mod
use fhash_module__ints_double
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                                    :: numNodes, targetNumGraftedChains
integer, intent(in), dimension(targetNumGraftedChains) :: gpid
integer                                                :: ii, jj

type(fhash_type__ints_double), intent(inout) :: elemcon

real(8), intent(in), dimension(targetNumGraftedChains)            :: gp_init_value
real(8), intent(in), dimension(numNodes)                          :: ww
real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qmx_interp_mg
real(8), intent(in), dimension(numEdwPointsGrafted+1)             :: ds_gr_ed, xs_gr_ed
real(8), intent(in), dimension(numConvolPointsGrafted+1)          :: xs_gr_conv, coeff_gr_conv
real(8), intent(out), dimension(numNodes,targetNumGraftedChains)  :: phi_gr_indiv
real(8), dimension(2,numNodes)                                    :: qgr
real(8), dimension(numEdwPointsGrafted+1,numNodes)                :: qgr_final
real(8), dimension(numConvolPointsGrafted+1,numNodes)             :: qgr_interp

integer:: gpIndex
!------------------------------------------------------------------------------------------------------!
write(6,'(2X,A43)')adjl("Computing indiv profiles of grafted chains.",43)

call FemMatrixAssemble(rg2OfGraftedMonomer, ww)

if (exportAllGraftedChains.eq.1) then ! export profile of all grafted chains separately
  do ii = 1, targetNumGraftedChains
    qgr       = 0.0d0
    qgr_final = 0.0d0

    qgr(1,gpid(ii))       = gp_init_value(ii)
    qgr_final(1,gpid(ii)) = gp_init_value(ii)

    write(6, '(2X,A21,1X,I3,1X,A8,1X,I7,1X,A2,1X)', advance='no') "Grafting point index:", ii, "with id:", gpid(ii), "->"
    call SolverEdwards(ds_gr_ed, numEdwPointsGrafted, qgr, qgr_final, nodeBelongsToDirichletFace, elemcon)

    do jj = 1, numNodes
      call interp_linear(1, numEdwPointsGrafted+1, xs_gr_ed, qgr_final(:,jj), numConvolPointsGrafted+1, xs_gr_conv, qgr_interp(:,jj))
    enddo

    call ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeff_gr_conv, qgr_interp, qmx_interp_mg, phi_gr_indiv(:,ii))
  enddo
elseif (exportAllGraftedChains.eq.0) then ! export profile of specific grafted chains separately
  do ii = 1, numGraftedChainsToExport
    qgr       = 0.0d0
    qgr_final = 0.0d0

    gpIndex = gpIndexToExport(ii)

    qgr(1,gpid(gpIndex))       = gp_init_value(gpIndex)
    qgr_final(1,gpid(gpIndex)) = gp_init_value(gpIndex)

    write(6, '(2X,A21,1X,I3,1X,A8,1X,I7,1X,A2,1X)', advance='no') "Grafting point index:", gpIndex, "with id:", gpid(gpIndex), "->"
    call SolverEdwards(ds_gr_ed, numEdwPointsGrafted, qgr, qgr_final, nodeBelongsToDirichletFace, elemcon)

    do jj = 1, numNodes
      call interp_linear(1, numEdwPointsGrafted+1, xs_gr_ed, qgr_final(:,jj), numConvolPointsGrafted+1, xs_gr_conv, qgr_interp(:,jj))
    enddo

    call ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeff_gr_conv, qgr_interp, qmx_interp_mg, phi_gr_indiv(:,gpIndex))
  enddo
elseif (exportAllGraftedChains.eq.2) then ! export profile of specific grafted chains as a sum
  qgr       = 0.0d0
  qgr_final = 0.0d0

  do ii = 1, numGraftedChainsToExport
    gpIndex = gpIndexToExport(ii)

    qgr(1,gpid(gpIndex))       = gp_init_value(gpIndex)
    qgr_final(1,gpid(gpIndex)) = gp_init_value(gpIndex)
  enddo

  write(6, '(2X,A40)') "Some of specifiec grafted chains:"
  call SolverEdwards(ds_gr_ed, numEdwPointsGrafted, qgr, qgr_final, nodeBelongsToDirichletFace, elemcon)

  do jj = 1, numNodes
    call interp_linear(1, numEdwPointsGrafted+1, xs_gr_ed, qgr_final(:,jj), numConvolPointsGrafted+1, xs_gr_conv, qgr_interp(:,jj))
  enddo

  call ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeff_gr_conv, qgr_interp, qmx_interp_mg, phi_gr_indiv(:,1))
endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine ComputeIndivProfile
