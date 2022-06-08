!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportEndMiddleProfile(ns, q1_final, qmx_final, chain_type, numNodes, nodeCoord)
!----------------------------------------------------------------------------------------------------------!
use write_helper_mod, only : adjl
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns, numNodes
integer             :: kk, jj, ns_middle

real(8), intent(in), dimension(ns,numNodes) :: q1_final, qmx_final
real(8), intent(in), dimension(3,numNodes)  :: nodeCoord
real(8), dimension(numNodes)                :: phi_middle, phi_end

character(2), intent(in) :: chain_type
character(80)            :: filename = ''
!----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting end/middle density profiles.",40)

phi_middle = 0.0d0
phi_end    = 0.0d0

do kk = 1, numNodes
  phi_end(kk) = (q1_final(ns,kk) * qmx_final(1,kk))
enddo

ns_middle = ns/2

do kk = 1, numNodes
  phi_middle(kk) = q1_final(ns_middle,kk) * qmx_final(ns_middle+1,kk)
enddo

write(filename,'("o.phi_end_middle_",A2)') chain_type

open(unit=122, file=filename)
write(122,'(5(2X,A16))') 'x', 'y', 'z', "phi_end", "phi_middle"
do kk = 1, numNodes
  write (122,'(5(2X,E16.9))') (nodeCoord(jj,kk), jj=1,3), phi_end(kk), phi_middle(kk)
enddo
close(unit=122)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine ExportEndMiddleProfile
