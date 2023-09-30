!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportFieldAscii(ww, wwNew, wwMix)
!-----------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: numNodes, nodeCoord
use iofiles_mod,      only: IO_field
use write_helper_mod, only: adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: node

real(8), intent(in), dimension(numNodes) :: ww, wwNew, wwMix
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting field.",40)
open (unit=120, file = IO_field)
write(120,'(7(A19))') "np", "x", "y", "z", "ww", "ww_new", "ww_mix"
do node = 1, numNodes
  write(120,'(I19,6(E19.9E3))') node, nodeCoord(1,node), nodeCoord(2,node), nodeCoord(3,node), ww(node), wwNew(node), wwMix(node)
enddo
close(120)
!-----------------------------------------------------------------------------------------------------------!
end subroutine ExportFieldAscii
