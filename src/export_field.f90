!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_field(ww, ww_new, ww_mix)
!-----------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: numNodes, nodeCoord
use iofiles_mod,      only: field_profile
use write_helper_mod, only: adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: kk

real(8), intent(in), dimension(numNodes) :: ww, ww_new, ww_mix
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting field.",40)
open (unit=120, file = field_profile)
write(120,'(7(A19))') "np", "x", "y", "z", "ww", "ww_new", "ww_mix"
do kk = 1, numNodes
  write(120,'(I19,6(E19.9E3))') kk, nodeCoord(1,kk), nodeCoord(2,kk), nodeCoord(3,kk), ww(kk), ww_new(kk), ww_mix(kk)
enddo
close(120)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_field
