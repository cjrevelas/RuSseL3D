!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module hist_mod
!----------------------------------------------------------------------------------------------------------------------------!
integer                                :: numBins
integer, allocatable, dimension(:,:,:) :: planarCellId
integer, allocatable, dimension(:,:)   :: sphericalCellId

real(8)                                :: binLength
real(8), allocatable, dimension(:,:,:) :: distanceFromFace, planarCellVolume
real(8), allocatable, dimension(:,:)   :: distanceFromNanop, sphericalCellVolume
!----------------------------------------------------------------------------------------------------------------------------!
end module hist_mod
