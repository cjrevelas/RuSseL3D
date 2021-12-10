!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module kcw
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer                              :: NNZ
real(8), allocatable, dimension(:)   :: rdiag1

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

type(full_matrix) :: F_m
!--------------------------------------------------------------------!
end module kcw
