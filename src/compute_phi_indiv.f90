subroutine compute_phi_indiv(numnp, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, wa, &
     &                num_gpoints, gpid, gp_init_value, phia_gr_indiv)
!------------------------------------------------------------------------------------------------------!
use geometry,     only: node_in_q0_face
use write_helper, only: adjl
use parser_vars,  only: ns_gr_conv, ns_gr_ed, chainlen_gr, mumps_matrix_type, rg2_per_mon_gr
use error_handing
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: numnp, num_gpoints
integer, intent(in), dimension(num_gpoints) :: gpid
integer                                     :: ii, jj

real(8), intent(in), dimension(num_gpoints)        :: gp_init_value
real(8), intent(in), dimension(numnp)              :: wa
real(8), intent(in), dimension(ns_gr_conv+1,numnp) :: qmx_interp_mg
real(8), intent(in), dimension(ns_gr_ed+1)         :: ds_gr_ed, xs_gr_ed
real(8), intent(in), dimension(ns_gr_conv+1)       :: xs_gr_conv, coeff_gr_conv
real(8), intent(out), dimension(numnp,num_gpoints) :: phia_gr_indiv
real(8), dimension(2,numnp)                        :: qgr
real(8), dimension(ns_gr_ed+1,numnp)               :: qgr_final
real(8), dimension(ns_gr_conv+1,numnp)             :: qgr_interp
!------------------------------------------------------------------------------------------------------!
write(6,'(6X,A40)')adjl("*computing densities of grafted chains..",40)

call matrix_assemble(Rg2_per_mon_gr, wa)

do ii = 1, num_gpoints
    qgr       = 0.d0
    qgr_final = 0.d0

    qgr(1,gpid(ii))       = gp_init_value(ii)
    qgr_final(1,gpid(ii)) = gp_init_value(ii)

    write(6, '(9x,A9,I10)', advance='no') "..gp Id: ", gpid(ii)
    call edwards(ds_gr_ed, ns_gr_ed, mumps_matrix_type, qgr, qgr_final, node_in_q0_face)

    do jj = 1, numnp
        call interp_linear(1, ns_gr_ed+1, xs_gr_ed, qgr_final(:,jj), ns_gr_conv+1, xs_gr_conv, qgr_interp(:,jj))
    enddo

    call convolution(numnp, chainlen_gr, ns_gr_conv, coeff_gr_conv, qgr_interp, qmx_interp_mg, phia_gr_indiv(:,ii))
enddo

return
!------------------------------------------------------------------------------------------------------!
end subroutine compute_phi_indiv
