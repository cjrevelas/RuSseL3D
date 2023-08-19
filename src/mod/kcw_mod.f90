!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module kcw_mod
!----------------------------------------------------------------------------------------------------------------------------!
integer :: numNonZeroEntries

real(8), allocatable, dimension(:) :: rdiag1

type mumpsMatrix
  sequence
    double precision, dimension(:), pointer :: value
    integer, dimension(:), pointer          :: row, col
end type mumpsMatrix

type(mumpsMatrix) :: A_m

type fullMatrix
  sequence
    double precision, dimension(:), pointer :: g, rh, c, k, w
    integer, dimension(:), pointer          :: row, col
    logical, dimension(:), pointer          :: isZero
end type fullMatrix

type(fullMatrix) :: F_m
!----------------------------------------------------------------------------------------------------------------------------!
end module kcw_mod
