!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module delta_mod
!----------------------------------------------------------------------------------------------------------------------------!
integer                            :: targetNumGraftedChains
integer, allocatable, dimension(:) :: graftPointId

real(8), allocatable, dimension(:) :: deltaNumerical, graftPointValue
!----------------------------------------------------------------------------------------------------------------------------!
end module delta_mod
