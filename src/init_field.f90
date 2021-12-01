subroutine init_field(field_in_filename, Ufield, wa)
!------------------------------------------------------------------------------------------------------!
use parser_vars
use kcw
use error_handing
use write_helper
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: k1

character(20), intent(in) :: field_in_filename

real(8), intent(out), dimension(numnp) :: Ufield, wa
real(8)                                :: surf_pot, distance
!------------------------------------------------------------------------------------------------------!
open(unit=211, file = 'Usolid.out.txt')
do k1 = 1, numnp
    distance   = xc(1,k1)
    Ufield(k1) = surf_pot(distance)
    write(211,('(E16.9,2X,E19.9)')) distance, Ufield(k1)
    if (Ufield(k1).ne.Ufield(k1)) then
        Ufield(k1) = 0.d0
        write(ERROR_MESSAGE,'(''Hamaker assumed a NaN value for x = '',E16.9,''. NaN was changed to '',E16.9)') distance, Ufield(k1)
        call exit_with_error(0,2,0,ERROR_MESSAGE)
    endif
enddo
close(211)
!*******************************************************************!
!                             READ FIELD                            !
!*******************************************************************!
if (readfield.eq.1) then
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
else
    if (scheme_type.eq.2) then
        do k1 = 1, numnp
            if (elem_in_q0_face(k1)) then
                wa(k1) = -kapa
            else
                wa(k1) = 0.d0
            endif
        enddo
    else
        wa = 0.d0
    endif
endif
!------------------------------------------------------------------------------------------------------!
end subroutine init_field
