subroutine export_field_bin(wa, numnp, iter)
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in)                   :: numnp, iter

real(8), intent(in), dimension(numnp) :: wa

character(40) :: field_filename_aux = ''
!-----------------------------------------------------------------------------------------------------------!
write(field_filename_aux,'("field",I4.4,".out.bin")') iter
open(unit=655, file = field_filename_aux, form='unformatted')
write(655) wa
close(655)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_field_bin
