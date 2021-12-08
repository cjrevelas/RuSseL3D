subroutine find_delta(numnp, qm_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, koeff_gr_conv, wa, &
     &                num_gpoints, gpid, delta_numer, gp_init_value)
!------------------------------------------------------------------------------------------------------!
use error_handing
use iofiles
use parser_vars
use constants
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: numnp, num_gpoints
integer, intent(in), dimension(num_gpoints) :: gpid
integer                                     :: i1, j1, iog

real(8), dimension(num_gpoints), intent(out) :: delta_numer
real(8), dimension(num_gpoints), intent(in)  :: gp_init_value

real(8), allocatable, dimension(:)                 :: delta_anal
real(8), dimension(numnp), intent(in)              :: wa
real(8), dimension(ns_gr_conv+1,numnp), intent(in) :: qm_interp_mg
real(8), dimension(ns_gr_ed+1), intent(in)         :: ds_gr_ed, xs_gr_ed
real(8), dimension(ns_gr_conv+1), intent(in)       :: xs_gr_conv, koeff_gr_conv
real(8), dimension(numnp)                          :: phia_gr
real(8), dimension(2,numnp)                        :: qgr
real(8), dimension(ns_gr_ed+1,numnp)               :: qgr_final
real(8), dimension(ns_gr_conv+1,numnp)             :: qgr_interp
real(8)                                            :: aux1 = 0.d0, aux2 = 0.d0, nch_gr = 0.d0, node_vol = 0.d0
real(8)                                            :: initValue = 0.d0
!------------------------------------------------------------------------------------------------------!
allocate(delta_anal(num_gpoints))

delta_anal = 0.d0

!analytic delta calculation
do i1 = 1, num_gpoints
    phia_gr           = 0.d0
    phia_gr(gpid(i1)) = 1.d0
    call spat_3d(phia_gr, node_vol, aux1, aux2)
    delta_anal(i1) = 1.d0 / node_vol * 1.e+30
enddo

!numerical delta calculation (simulation)
call matrix_assemble(Rg2_per_mon_gr, wa)

do i1 = 1, num_gpoints
    qgr       = 0.d0
    qgr_final = 0.d0

    initValue = delta_anal(i1) * chainlen_gr &
                               * 1.d0 / (qm_interp_mg(ns_gr_conv+1,gpid(i1)) * (rho_mol_bulk * n_avog))

    qgr(1,gpid(i1))       = initValue
    qgr_final(1,gpid(i1)) = initValue

    write(6, '(A30,I10)', advance='no') "Computing delta for node: ", gpid(i1)
    call edwards(ds_gr_ed, ns_gr_ed, mumps_matrix_type, qgr, qgr_final)

    do j1 = 1, numnp
        call interp_linear(1, ns_gr_ed+1, xs_gr_ed, qgr_final(:,j1), ns_gr_conv+1, xs_gr_conv, qgr_interp(:,j1))
    enddo

    call convolution(numnp, chainlen_gr, ns_gr_conv, koeff_gr_conv, qgr_interp, qm_interp_mg, phia_gr)

    call get_nchains(numnp, chainlen_gr, rho_mol_bulk, phia_gr, nch_gr)

    delta_numer(i1) = delta_anal(i1) / nch_gr
enddo

!update gnodes.in.lammpstrj input file
iog = 19
open(unit=iog, file=gp_filename)
write(iog,'(A14)') "ITEM: TIMESTEP"
write(iog,'(A1)')  '0'
write(iog,'(A21)') "ITEM: NUMBER OF ATOMS"
write(iog,'(I10)')  num_gpoints
write(iog,'(A25)') "ITEM: BOX BOUNDS PP PP PP"
do i1 = 1, 3
    write(iog,'(2(F19.15))') box_lo(i1), box_hi(i1)
enddo
write(iog,'(A28)') "ITEM: ATOMS id type xu yu zu"

do i1 = 1, num_gpoints
    write(iog,'(I5,2(1X,E16.9))') gpid(i1), gp_init_value(i1), delta_numer(i1)
enddo
close(iog)

open(unit=5, file = delta_out)
write(5,'(A5,A16,A20,A17)') "gpid", "qm(rgi,Ng)", "Numerical delta", "Analytic delta"
do i1 = 1, num_gpoints
    write(5,'(I5,3(1X,E18.9))') gpid(i1), qm_interp_mg(ns_gr_conv+1,gpid(i1)), delta_numer(i1), delta_anal(i1)
enddo
close(5)

deallocate(delta_anal)
!------------------------------------------------------------------------------------------------------!
end subroutine find_delta
