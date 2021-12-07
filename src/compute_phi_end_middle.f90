subroutine compute_phi_end_middle(ns, q1_final, qmatrix_final, chain_type)
!----------------------------------------------------------------------------------------------------------!
use constants
use geometry
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: kk, jj, ns_middle

real(8), intent(in), dimension(numnp,ns) :: q1_final, qmatrix_final
real(8), dimension(numnp)                :: phi_middle, phi_end

character, intent(in) :: chain_type
character(40)         :: filename = ''
!----------------------------------------------------------------------------------------------------------!
phi_middle = 0.d0
phi_end    = 0.d0

do kk = 1, numnp
    phi_end(kk) = (q1_final(kk,ns) * qmatrix_final(kk,1))
enddo

ns_middle = ns/2

do kk = 1, numnp
    phi_middle(kk) = q1_final(kk,ns_middle) * qmatrix_final(kk,ns_middle+1)
enddo

write(filename,'("phi_end_middle_",A6,".out.txt")') chain_type

open(unit=122, file=filename)
write(122,'(5(2X,A16))') 'x', 'y', 'z', "phi_end", "phi_middle"
do kk = 1, numnp
    write (122,'(5(2X,E16.9))') (xc(jj,kk), jj=1,3), phi_end(kk), phi_middle(kk)
enddo
close(unit=122)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_phi_end_middle
