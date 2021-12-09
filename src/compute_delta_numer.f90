subroutine compute_delta_numer(numnp, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, wa, &
     &                num_gpoints, gpid, delta_numer, volnp)
!------------------------------------------------------------------------------------------------------!
use geometry,     only: node_in_q0_face
use parser_vars,  only: ns_gr_conv, ns_gr_ed, chainlen_gr, mumps_matrix_type, rg2_per_mon_gr, rho_mol_bulk
use constants,    only: n_avog, m3_to_A3
use write_helper, only: adjl
use error_handing
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: numnp, num_gpoints
integer, intent(in), dimension(num_gpoints) :: gpid
integer                                     :: ii, jj

real(8), intent(in), dimension(numnp)              :: wa, volnp
real(8), intent(in), dimension(ns_gr_conv+1,numnp) :: qmx_interp_mg
real(8), intent(in), dimension(ns_gr_ed+1)         :: ds_gr_ed, xs_gr_ed
real(8), intent(in), dimension(ns_gr_conv+1)       :: xs_gr_conv, coeff_gr_conv
real(8), intent(out), dimension(num_gpoints)       :: delta_numer
real(8), dimension(num_gpoints)                    :: delta_anal
real(8), dimension(numnp)                          :: phia_gr
real(8), dimension(2,numnp)                        :: qgr
real(8), dimension(ns_gr_ed+1,numnp)               :: qgr_final
real(8), dimension(ns_gr_conv+1,numnp)             :: qgr_interp
real(8)                                            :: nch_gr = 0.d0, initValue = 0.d0
!------------------------------------------------------------------------------------------------------!
write(6,'(6X,A40)')adjl("..updating delta of grafted chains:",40)

delta_anal = 0.d0

!analytic delta calculation
do ii = 1, num_gpoints
    delta_anal(ii) = 1.d0 / volnp(gpid(ii)) * m3_to_A3
enddo

!numerical delta calculation (simulation)
call matrix_assemble(Rg2_per_mon_gr, wa)

do ii = 1, num_gpoints
    qgr       = 0.d0
    qgr_final = 0.d0

    initValue = delta_anal(ii) * chainlen_gr &
                               * 1.d0 / (qmx_interp_mg(ns_gr_conv+1,gpid(ii)) * (rho_mol_bulk * n_avog))

    qgr(1,gpid(ii))       = initValue
    qgr_final(1,gpid(ii)) = initValue

    write(6, '(9x,A9,I10)', advance='no') "..gp id: ", gpid(ii)
    call edwards(ds_gr_ed, ns_gr_ed, mumps_matrix_type, qgr, qgr_final, node_in_q0_face)

    do jj = 1, numnp
        call interp_linear(1, ns_gr_ed+1, xs_gr_ed, qgr_final(:,jj), ns_gr_conv+1, xs_gr_conv, qgr_interp(:,jj))
    enddo

    call convolution(numnp, chainlen_gr, ns_gr_conv, coeff_gr_conv, qgr_interp, qmx_interp_mg, phia_gr)

    call get_nchains(numnp, chainlen_gr, rho_mol_bulk, phia_gr, nch_gr)

    delta_numer(ii) = delta_anal(ii) / nch_gr
enddo
!------------------------------------------------------------------------------------------------------!
end subroutine compute_delta_numer