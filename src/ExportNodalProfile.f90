!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportNodalProfile(phiMatrix, phiGrafted, numNodes, nodeCoord, nodeVolume)
!-----------------------------------------------------------------------------------------------------------!
use iofiles_mod,      only: IO_nodalProfile
use write_helper_mod, only: adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numNodes
integer             :: node

real(8), intent(in), dimension(numNodes)   :: nodeVolume
real(8), intent(in), dimension(numNodes)   :: phiMatrix, phiGrafted
real(8), intent(in), dimension(3,numNodes) :: nodeCoord
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("****************************************",40)
write(6,'(2X,A40)')adjl("Exporting density profiles.",40)
open (unit=120, file = IO_nodalProfile)
write(120,'(7(A19))') "np", "x", "y", "z", "phi_m", "phi_g", "vol"
do node = 1, numNodes
  write(120,'(I19,6(E19.9E3))') node, nodeCoord(1,node), nodeCoord(2,node), nodeCoord(3,node), phiMatrix(node), phiGrafted(node), nodeVolume(node)
enddo
close(120)
!-----------------------------------------------------------------------------------------------------------!
end subroutine ExportNodalProfile
