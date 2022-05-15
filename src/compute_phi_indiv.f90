!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_phi_indiv(numnp, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, ww, &
     &                targetNumGraftedChains, gpid, gp_init_value, phi_gr_indiv)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: nodeBelongsToDirichletFace
use write_helper_mod, only: adjl
use parser_vars_mod,  only: ns_gr_conv, ns_gr_ed, lengthGrafted, mumpsMatrixType, rg2OfGraftedMonomer
use error_handing_mod
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: numnp, targetNumGraftedChains
integer, intent(in), dimension(targetNumGraftedChains) :: gpid
integer                                     :: ii, jj

real(8), intent(in), dimension(targetNumGraftedChains)        :: gp_init_value
real(8), intent(in), dimension(numnp)              :: ww
real(8), intent(in), dimension(ns_gr_conv+1,numnp) :: qmx_interp_mg
real(8), intent(in), dimension(ns_gr_ed+1)         :: ds_gr_ed, xs_gr_ed
real(8), intent(in), dimension(ns_gr_conv+1)       :: xs_gr_conv, coeff_gr_conv
real(8), intent(out), dimension(numnp,targetNumGraftedChains) :: phi_gr_indiv
real(8), dimension(2,numnp)                        :: qgr
real(8), dimension(ns_gr_ed+1,numnp)               :: qgr_final
real(8), dimension(ns_gr_conv+1,numnp)             :: qgr_interp
!------------------------------------------------------------------------------------------------------!
write(6,'(2X,A43)')adjl("Computing indiv profiles of grafted chains.",43)

call fem_matrix_assemble(rg2OfGraftedMonomer, ww)

do ii = 1, targetNumGraftedChains
    qgr       = 0.0d0
    qgr_final = 0.0d0

    qgr(1,gpid(ii))       = gp_init_value(ii)
    qgr_final(1,gpid(ii)) = gp_init_value(ii)

    write(6, '(2x,A19,I7,A3)', advance='no') "Grafting point id: ", gpid(ii), " ->"
    call solver_edwards(ds_gr_ed, ns_gr_ed, mumpsMatrixType, qgr, qgr_final, nodeBelongsToDirichletFace)

    do jj = 1, numnp
        call interp_linear(1, ns_gr_ed+1, xs_gr_ed, qgr_final(:,jj), ns_gr_conv+1, xs_gr_conv, qgr_interp(:,jj))
    enddo

    call contour_convolution(numnp, lengthGrafted, ns_gr_conv, coeff_gr_conv, qgr_interp, qmx_interp_mg, phi_gr_indiv(:,ii))
enddo

return
!------------------------------------------------------------------------------------------------------!
end subroutine compute_phi_indiv
