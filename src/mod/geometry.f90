module geometry
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer                              :: all_el, nel, ndm, numnp, numel
integer, allocatable, dimension(:)   :: con_l2
integer, allocatable, dimension(:,:) :: ix

logical, dimension(3,2)              :: is_dir_face
logical, allocatable, dimension(:)   :: elem_in_q0_face

real(8), allocatable, dimension(:,:) :: xc
!--------------------------------------------------------------------!
end module geometry
