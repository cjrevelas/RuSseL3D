subroutine export_phi_nodal(phia_mx, phia_gr, numnp, xc, volnp)
!-----------------------------------------------------------------------------------------------------------!
use iofiles, only: phi_nodal
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numnp
integer             :: kk

real(8), intent(in), dimension(numnp)   :: volnp
real(8), intent(in), dimension(numnp)   :: phia_mx, phia_gr
real(8), intent(in), dimension(3,numnp) :: xc
!-----------------------------------------------------------------------------------------------------------!
open (unit=120, file = phi_nodal)
write(120,'(7(A19))') "np", "x", "y", "z", "phi_m", "phi_g", "vol"
do kk = 1, numnp
   write(120,'(I19,6(E19.9E3))') kk, xc(1,kk), xc(2,kk), xc(3,kk), phia_mx(kk), phia_gr(kk), volnp(kk)
enddo
close(120)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_phi_nodal
