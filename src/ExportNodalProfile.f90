!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportNodalProfile(phi_mx, phi_gr, numNodes, nodeCoord, volnp)
!-----------------------------------------------------------------------------------------------------------!
use iofiles_mod,      only: IO_nodalProfile
use write_helper_mod, only: adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numNodes
integer             :: kk

real(8), intent(in), dimension(numNodes)   :: volnp
real(8), intent(in), dimension(numNodes)   :: phi_mx, phi_gr
real(8), intent(in), dimension(3,numNodes) :: nodeCoord
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("****************************************",40)
write(6,'(2X,A40)')adjl("Exporting density profiles.",40)
open (unit=120, file = IO_nodalProfile)
write(120,'(7(A19))') "np", "x", "y", "z", "phi_m", "phi_g", "vol"
do kk = 1, numNodes
  write(120,'(I19,6(E19.9E3))') kk, nodeCoord(1,kk), nodeCoord(2,kk), nodeCoord(3,kk), phi_mx(kk), phi_gr(kk), volnp(kk)
enddo
close(120)
!-----------------------------------------------------------------------------------------------------------!
end subroutine ExportNodalProfile
