module kcw
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: non_zero
integer :: all_el
integer :: min_el

integer, allocatable, dimension(:) :: con_l2
integer, allocatable, dimension(:,:) :: connectivity, con_l
double precision, allocatable :: k(:,:), c(:,:), w(:,:)
double precision, allocatable :: g(:,:), rdiag1(:), rh(:,:)

type mum_matrix
    sequence
    double precision, dimension(:), pointer :: value
    integer, dimension(:), pointer :: row, col
    !integer :: col
    !integer :: row
    !double precision :: m
end type

type(mum_matrix) :: g_m
type(mum_matrix) :: rh_m
type(mum_matrix) :: c_m
type(mum_matrix) :: c_min
type(mum_matrix) :: A_m
type(mum_matrix) :: k_m
type(mum_matrix) :: k_min
type(mum_matrix) :: w_m
type(mum_matrix) :: w_min
type(mum_matrix) :: g_min
type(mum_matrix) :: rh_min
!--------------------------------------------------------------------!        
end module kcw 