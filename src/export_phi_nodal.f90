!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_phi_nodal(phi_mx, phi_gr, numnp, xc, volnp)
!-----------------------------------------------------------------------------------------------------------!
use iofiles_mod,      only: phi_nodal
use write_helper_mod, only: adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numnp
integer             :: kk

real(8), intent(in), dimension(numnp)   :: volnp
real(8), intent(in), dimension(numnp)   :: phi_mx, phi_gr
real(8), intent(in), dimension(3,numnp) :: xc
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("****************************************",40)
write(6,'(2X,A40)')adjl("Exporting density profiles.",40)
open (unit=120, file = phi_nodal)
write(120,'(7(A19))') "np", "x", "y", "z", "phi_m", "phi_g", "vol"
do kk = 1, numnp
   write(120,'(I19,6(E19.9E3))') kk, xc(1,kk), xc(2,kk), xc(3,kk), phi_mx(kk), phi_gr(kk), volnp(kk)
enddo
close(120)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_phi_nodal
