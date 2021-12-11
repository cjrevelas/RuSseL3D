!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine get_delta_numer(numnp, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, wa, num_gpoints, gpid, delta_numer, volnp)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: node_belongs_to_dirichlet_face
use parser_vars_mod,  only: ns_gr_conv, ns_gr_ed, chainlen_gr, mumps_matrix_type, rg2_per_mon_gr, rho_mol_bulk
use constants_mod,    only: n_avog, m3_to_A3
use write_helper_mod, only: adjl
use error_handing_mod
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
write(6,'(2X,A40)')adjl("****************************************",40)
write(6,'(2X,A40)')adjl("Updating delta of grafted chains:",40)

delta_anal = 0.d0

! Analytic delta calculation
do ii = 1, num_gpoints
    delta_anal(ii) = 1.d0 / volnp(gpid(ii)) * m3_to_A3
enddo

! Numerical delta calculation (i.e., through solution of the Edwards equation)
call fem_matrix_assemble(Rg2_per_mon_gr, wa)

do ii = 1, num_gpoints
    qgr       = 0.d0
    qgr_final = 0.d0

    initValue = delta_anal(ii) * chainlen_gr * 1.d0 / (qmx_interp_mg(ns_gr_conv+1,gpid(ii)) * (rho_mol_bulk * n_avog))

    qgr(1,gpid(ii))       = initValue
    qgr_final(1,gpid(ii)) = initValue

    write(6, '(2X,A19,I7,A3)', advance='no') "Grafting point id: ", gpid(ii), " ->"
    call solver_edwards(ds_gr_ed, ns_gr_ed, mumps_matrix_type, qgr, qgr_final, node_belongs_to_dirichlet_face)

    do jj = 1, numnp
        call interp_linear(1, ns_gr_ed+1, xs_gr_ed, qgr_final(:,jj), ns_gr_conv+1, xs_gr_conv, qgr_interp(:,jj))
    enddo

    call contour_convolution(numnp, chainlen_gr, ns_gr_conv, coeff_gr_conv, qgr_interp, qmx_interp_mg, phia_gr)

    call compute_number_of_chains(numnp, chainlen_gr, rho_mol_bulk, phia_gr, nch_gr)

    delta_numer(ii) = delta_anal(ii) / nch_gr
enddo
write(6,'(2X,A40)')adjl("****************************************",40)
!------------------------------------------------------------------------------------------------------!
end subroutine get_delta_numer
