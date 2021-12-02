subroutine init_field(field_in_filename, Ufield, wa)
!------------------------------------------------------------------------------------------------------!
use parser_vars
use geometry
use error_handing
use write_helper
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: k1, m1, n1

character(20), intent(in) :: field_in_filename

real(8), intent(out), dimension(numnp) :: Ufield, wa
real(8)                                :: distance, r12, Urep, Uatt, Utot
!------------------------------------------------------------------------------------------------------!
open(unit=211, file = 'Usolid.out.txt')
Ufield = 0.d0
do k1 = 1, numnp
    ! loop over all dirichlet faces
    do m1 = 1, 3
        do n1 = 1, 2
            if (is_dir_face(m1,n1)) then

                if (n1.eq.1) then
                    distance = xc(m1,k1) - box_lo(m1)
                elseif (n1.eq.2) then
                    distance = box_hi(m1) - xc(m1,k1)
                endif

                call surf_pot(Temp, distance, sphere_radius, rho_0, sigma1, sigma2, Aps, Asio2, r12, Urep, Uatt, Utot)
                Ufield(k1) = Ufield(k1) + Utot

                if (Ufield(k1).ne.Ufield(k1)) then
                    write(ERROR_MESSAGE,'(''Hamaker assumed a NaN value for x = '',E16.9,''. &
                                        & NaN was changed to '',E16.9)') distance, Ufield(k1)
                    call exit_with_error(1,2,1,ERROR_MESSAGE)
                endif

                if (distance.lt.0.d0) then
                    write(ERROR_MESSAGE,'(''Hamaker distance smaller than zero! ('',E16.9,'').'')') distance
                    call exit_with_error(1,2,1,ERROR_MESSAGE)
                endif

                write(211,'(E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3)') distance, r12, Uatt, Urep, Utot
            endif
        enddo
    enddo
    !Nanoparticle section
    !       distance = dsqrt(xc(1,k1)**2 + xc(2,k1)**2 + xc(3,k1)**2)
enddo

close(211)
!*******************************************************************!
!                             READ FIELD                            !
!*******************************************************************!
if (field_init_scheme.eq.0) then
    wa = 0.d0
elseif (field_init_scheme.eq.1) then
    write(iow,'(/A40,5x,A12)')adjl('*Reading field from file:',40),field_in_filename
    write(6  ,'(/A40,5x,A12)')adjl('*Reading field from file:',40),field_in_filename

    INQUIRE(FILE=field_in_filename, EXIST=FILE_EXISTS)

    if (FILE_EXISTS) then
        open(unit=655, file = field_in_filename, Form='unformatted')
    else
        write(ERROR_MESSAGE,'(''File '',A15,'' does not exist!'')')field_in_filename
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    read(655) wa
    close(655)
elseif (field_init_scheme.eq.1) then
    do k1 = 1, numnp
        if (elem_in_q0_face(k1)) then
            wa(k1) = -kapa
        else
            wa(k1) = 0.d0
        endif
    enddo
endif
!------------------------------------------------------------------------------------------------------!
end subroutine init_field
