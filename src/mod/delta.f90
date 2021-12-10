!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module delta
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer                            :: num_gpoints
integer, allocatable, dimension(:) :: gpid
real(8), allocatable, dimension(:) :: delta_numer, gp_init_value
!-----------------------------------------------------------------------------------------------------------!
end module delta
