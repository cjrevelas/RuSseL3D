!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_field_bin(ww, numNodes, iter)
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numNodes, iter

real(8), intent(in), dimension(numNodes) :: ww

character(40) :: field_filename_aux = ''
!-----------------------------------------------------------------------------------------------------------!
write(field_filename_aux,'("o.field",I4.4,".bin")') iter
open(unit=655, file = field_filename_aux, form='unformatted')
write(655) ww
close(655)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_field_bin
