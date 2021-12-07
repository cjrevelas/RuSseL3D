subroutine compute_phi_end_middle(ns, nx, rx, q1_final, qmatrix_final, chain_type)
!----------------------------------------------------------------------------------------------------------!
use constants
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns, nx
integer             :: kk, ns_middle

real(8), intent(in), dimension(0:nx)      :: rx
real(8), intent(in), dimension(0:nx,0:ns) :: q1_final, qmatrix_final
real(8), dimension(0:nx)                  :: phi_middle, phi_end

character(6), intent(in) :: chain_type
character(40)            :: filename = ''
!----------------------------------------------------------------------------------------------------------!
phi_middle = 0.d0
phi_end    = 0.d0

!start profile post-processing
do kk = 0, nx
    phi_end(kk) = (q1_final(kk,ns) * qmatrix_final(kk,0))
enddo

ns_middle = ns/2

do kk = 0, nx
    phi_middle(kk) = q1_final(kk,ns_middle) * qmatrix_final(kk,ns_middle +1)
enddo

write(filename,'(''o.phi_end_middle_'',A6,''.out.txt'')')chain_type

open(unit=122, file=filename)
write(122,'(3(2X,A16))') 'r', 'phi_end(r)', 'phi_middle(r)'
do kk = 0, nx
    write (122,'(3(2X,E16.9))') rx(kk), phi_end(kk), phi_middle(kk)
enddo
close(unit=122)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_phi_end_middle
