!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_field_bin(wa, numnp, iter)
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in)                   :: numnp, iter

real(8), intent(in), dimension(numnp) :: wa

character(40) :: field_filename_aux = ''
!-----------------------------------------------------------------------------------------------------------!
write(field_filename_aux,'("o.field",I4.4,".bin")') iter
open(unit=655, file = field_filename_aux, form='unformatted')
write(655) wa
close(655)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_field_bin
