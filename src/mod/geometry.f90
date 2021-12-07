module geometry
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer                              :: all_el, nel, ndm, numnp, numel, max_el_node
integer, allocatable, dimension(:)   :: con_l2, n_el_node
integer, allocatable, dimension(:,:) :: ix, el_node

logical, dimension(3,2)              :: is_dir_face
logical, allocatable, dimension(:)   :: elem_in_q0_face

real(8), allocatable, dimension(:,:) :: xc

real(8), parameter :: dx = 1.D-01
real(8), parameter :: dy = 1.D-01
real(8), parameter :: dz = 1.D-01
!--------------------------------------------------------------------!
end module geometry
