!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module kcw_mod
!----------------------------------------------------------------------------------------------------------------------------!
integer :: NNZ

real(8), allocatable, dimension(:) :: rdiag1

type mumps_matrix
  sequence
    double precision, dimension(:), pointer :: value
    integer, dimension(:), pointer          :: row, col
end type mumps_matrix

type(mumps_matrix) :: A_m

type full_matrix
  sequence
    double precision, dimension(:), pointer :: g, rh, c, k, w
    integer, dimension(:), pointer          :: row, col
    logical, dimension(:), pointer          :: is_zero
end type full_matrix

type(full_matrix) :: F_m
!----------------------------------------------------------------------------------------------------------------------------!
end module kcw_mod
