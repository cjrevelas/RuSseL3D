module geometry
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer                              :: all_el, nel, ndm, numnp, numel, max_el_node
integer, allocatable, dimension(:)   :: con_l2, n_el_node
integer, allocatable, dimension(:,:) :: ix, el_node

logical, dimension(3,2)              :: is_dir_face
logical, allocatable, dimension(:)   :: node_in_q0_face

real(8), allocatable, dimension(:,:) :: xc

real(8), dimension(3) :: box_lo, box_hi, box_len

real(8), parameter :: dx = 1.d-01
real(8), parameter :: dy = 1.d-01
real(8), parameter :: dz = 1.d-01
!--------------------------------------------------------------------!
end module geometry
