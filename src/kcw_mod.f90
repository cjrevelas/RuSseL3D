module kcw
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer                              :: NNZ
integer                              :: all_el
integer, allocatable, dimension(:)   :: con_l2
integer, allocatable, dimension(:,:) :: connectivity, con_l

logical, allocatable, dimension(:) :: elem_in_q0_face

real(8), allocatable, dimension(:) :: rdiag1

type mumps_matrix
    sequence
    double precision, dimension(:), pointer :: value
    integer, dimension(:), pointer :: row, col
end type

type(mumps_matrix) :: A_m

type full_matrix
    sequence
    double precision, dimension(:), pointer :: g, rh, c, k, w
    integer, dimension(:), pointer :: row, col
end type

!*************************************MESH**************************************!
integer                              :: nel, ndm, numnp, numel
integer, allocatable, dimension(:,:) :: ix
real(8), allocatable, dimension(:,:) :: xc

type(full_matrix) :: F_m
!--------------------------------------------------------------------!
end module kcw
