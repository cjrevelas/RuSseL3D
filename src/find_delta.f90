subroutine find_delta(numnp, ds_gr_conv, koeff_gr_conv, wa, num_gpoints)
!------------------------------------------------------------------------------------------------------!
use error_handing
use iofiles
use parser_vars
use constants
use delta
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)  :: numnp
integer, intent(out) :: num_gpoints
integer              :: i1, iog

real(8), allocatable, dimension(:)           :: delta_anal
real(8), dimension(numnp), intent(in)        :: wa
real(8), dimension(ns_gr_conv+1), intent(in) :: ds_gr_conv, koeff_gr_conv
real(8), dimension(numnp)                    :: phia_gr, phia_gr_delta
real(8), dimension(numnp,2)                  :: qm_delta, qgr_delta
real(8), dimension(numnp,ns_gr_conv+1)       :: qm_final_delta, qgr_final_delta
real(8)                                      :: aux1 = 0.d0, aux2 = 0.d0, nch_gr = 0.d0, node_vol = 0.d0
!------------------------------------------------------------------------------------------------------!
iog = 19
inquire(file = gp_filename, exist = file_exists)

if (file_exists) then
    open(unit=iog, file = gp_filename)
else
    write(ERROR_MESSAGE,'("File ",A15," does not exist!")') gp_filename
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

read(iog,*)
read(iog,*)
read(iog,*)
read(iog,*) num_gpoints
read(iog,*)
read(iog,*)
read(iog,*)
read(iog,*)
read(iog,*)

allocate(gpid(num_gpoints), delta_numer(num_gpoints), gp_init_value(num_gpoints))

gpid = 0

delta_numer   = 0.d0
gp_init_value = 0.d0

!specify grafting points
do i1 = 1, num_gpoints

    read(iog,*) gpid(i1), gp_init_value(i1), delta_numer(i1)

    if (gpid(i1) > numnp) then
        write(ERROR_MESSAGE,'("ID of grafted chain (",I10,") is larger from numnp (",I10,")")') gpid(i1), numnp
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
enddo

if (grafted_ic_from_delta.eq.1) then

    !numerical delta calculation (simulation)
    call matrix_assemble(Rg2_per_mon_gr, wa)

    do i1 = 1, numnp
       qm_delta(i1,1)       = 1.d0
       qm_final_delta(i1,1) = 1.d0
    enddo

    call edwards(ds_gr_conv, ns_gr_conv, mumps_matrix_type, qm_delta, qm_final_delta)

    do i1 = 1, num_gpoints
        qgr_delta       = 0.d0
        qgr_final_delta = 0.d0

        qgr_delta(gpid(i1),1)       = 1.d0
        qgr_final_delta(gpid(i1),1) = 1.d0

        call edwards(ds_gr_conv, ns_gr_conv, mumps_matrix_type, qgr_delta, qgr_final_delta)

        call convolution(numnp, chainlen_gr, ns_gr_conv, koeff_gr_conv, qgr_final_delta, qm_final_delta, phia_gr)

        call grafted_chains(numnp, chainlen_gr, rho_0, phia_gr, nch_gr)

        delta_numer(i1) =(rho_0 * avogadro_constant) * qm_final_delta(gpid(i1), ns_gr_conv+1) * 1.d0 / (nch_gr * chainlen_gr)
    enddo

    !analytic delta calculation
    allocate(delta_anal(num_gpoints))

    delta_anal = 0.d0

    do i1 = 1, num_gpoints
        phia_gr_delta           = 0.d0
        phia_gr_delta(gpid(i1)) = 1.d0

        call spat_3d(phia_gr_delta, node_vol, aux1, aux2)

        delta_anal(i1) = 1.d0 / node_vol * 1.e+30
    enddo

    !update gnodes.in.lammpstrj input file
    rewind(iog)
    write(iog,'(A14)') "ITEM: TIMESTEP"
    write(iog,'(A1)')  '0'
    write(iog,'(A21)') "ITEM: NUMBER OF ATOMS"
    write(iog,'(A1)')  '1'
    write(iog,'(A25)') "ITEM: BOX BOUNDS PP PP PP"
    do i1 = 1, 3
        write(iog,'(2(F19.15))') box_lo(i1), box_hi(i1)
    enddo
    write(iog,'(A28)') "ITEM: ATOMS id type xu yu zu"

    do i1 = 1, num_gpoints
        write(iog,'(I5,2(E16.9))') gpid(i1), gp_init_value(i1), delta_numer(i1)
    enddo
    close(iog)

    open(unit=5, file = delta_out)
    write(5,'(A5,A16,A20,A17)') "gpid", "qm(rgi,Ng)", "Numerical delta", "Analytic delta"
    do i1 = 1, num_gpoints
        write(5,'(I5,3(E18.9))') gpid(i1), qm_final_delta(gpid(i1),ns_gr_conv+1), delta_numer(i1), delta_anal(i1)
    enddo
    close(5)

    deallocate(delta_anal)
endif
!------------------------------------------------------------------------------------------------------!
end subroutine find_delta
