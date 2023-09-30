!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportFieldBinary(ww, numNodes, iter)
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numNodes, iter

real(8), intent(in), dimension(numNodes) :: ww

character(40) :: fieldFilename = ''
!-----------------------------------------------------------------------------------------------------------!
write(fieldFilename,'("o.field",I4.4,".bin")') iter
open(unit=655, file = fieldFilename, form='unformatted')
write(655) ww
close(655)
!-----------------------------------------------------------------------------------------------------------!
end subroutine ExportFieldBinary
