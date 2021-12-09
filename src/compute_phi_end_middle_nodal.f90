subroutine compute_phi_end_middle_nodal(ns, q1_final, qmx_final, chain_type, numnp, xc)
!----------------------------------------------------------------------------------------------------------!
use write_helper, only : adjl
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns, numnp
integer             :: kk, jj, ns_middle

real(8), intent(in), dimension(ns,numnp) :: q1_final, qmx_final
real(8), intent(in), dimension(3,numnp)  :: xc
real(8), dimension(numnp)                :: phi_middle, phi_end

character, intent(in) :: chain_type
character(80)         :: filename = ''
!----------------------------------------------------------------------------------------------------------!
write(6,'(6X,A40)')adjl("*phi end/middle..",40)

phi_middle = 0.d0
phi_end    = 0.d0

do kk = 1, numnp
    phi_end(kk) = (q1_final(ns,kk) * qmx_final(1,kk))
enddo

ns_middle = ns/2

do kk = 1, numnp
    phi_middle(kk) = q1_final(ns_middle,kk) * qmx_final(ns_middle+1,kk)
enddo

write(filename,'("o.phi_end_middle_",A1,".out.txt")') chain_type

open(unit=122, file=filename)
write(122,'(A11,2A18,A21,A20)') 'x', 'y', 'z', "phi_end", "phi_middle"
do kk = 1, numnp
    write (122,'(5(2X,E16.9))') (xc(jj,kk), jj=1,3), phi_end(kk), phi_middle(kk)
enddo
close(unit=122)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_phi_end_middle_nodal