!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module geometry
!--------------------------------------------------------------------!
use parser_vars, only: radius_np_eff, n_nanopart_faces, wall_distance
use constants, only: pi
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer                              :: all_el, nel, ndm, numnp, numel, max_el_node
integer, allocatable, dimension(:)   :: con_l2, n_el_node
integer, allocatable, dimension(:,:) :: ix, el_node

logical, dimension(3,2)              :: is_dir_face
logical, allocatable, dimension(:)   :: node_in_q0_face

real(8), allocatable, dimension(:,:) :: xc
real(8), dimension(3)                :: box_lo, box_hi, box_len
real(8), parameter                   :: dx = 1.d-01
real(8), parameter                   :: dy = 1.d-01
real(8), parameter                   :: dz = 1.d-01
!real(8)                              :: interf_area
!--------------------------------------------------------------------!
contains
real(8) function interf_area()
   integer :: ii, mm, nn
   interf_area = 0.d0

   do mm = 1, 3
       do nn = 1, 2
           if (is_dir_face(mm,nn)) then
               interf_area = interf_area + DABS((box_hi(1)-box_lo(1))*(box_hi(2)-box_lo(2)))
           endif
       enddo
   enddo
   do ii = 1, n_nanopart_faces
       interf_area = interf_area + 4.d0*pi*(radius_np_eff(ii)-wall_distance)**2
       !OLD: interf_area = interf_area + 4.d0*pi*(radius_np_eff(ii)-4.d0)**2 * 4.d0
   enddo

return
end function interf_area
!--------------------------------------------------------------------!
end module geometry
