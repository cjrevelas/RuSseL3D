!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshProfile()
!----------------------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod, only: profileDimensions
use geometry_mod,    only: numNodes, boxLow, boxHigh, nodeCoord
use iofiles_mod,     only: IO_meshProfile
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
integer :: node, bin, numBins
integer, allocatable, dimension(:) :: profile

real(8) :: binSize
!----------------------------------------------------------------------------------------------------------------------------------!
binSize = 0.5d0
numBins = NINT((boxHigh(profileDimensions) - boxLow(profileDimensions)) / binSize) + 1

allocate(profile(numBins))

profile = 0

do node = 1, numNodes
  bin          = NINT((nodeCoord(profileDimensions,node) - boxLow(profileDimensions)) / binSize) + 1
  profile(bin) = profile(bin) + 1
enddo

open(77, file = IO_meshProfile)
do bin = 1, numBins
  write(77,*) bin, profile(bin)
enddo
close(77)

! Deallocate memory
deallocate(profile)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshProfile
