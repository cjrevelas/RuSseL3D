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
integer :: ii, ibin, nbin
integer, allocatable, dimension(:) :: prof_1D_node

real(8) :: prof_bin
!----------------------------------------------------------------------------------------------------------------------------------!
prof_bin = 0.5d0
nbin     = NINT((boxHigh(profileDimensions) - boxLow(profileDimensions)) / prof_bin) + 1

allocate(prof_1D_node(nbin))

prof_1D_node=0
do ii = 1, numNodes
  ibin               = NINT((nodeCoord(profileDimensions,ii) - boxLow(profileDimensions))/prof_bin) + 1
  prof_1D_node(ibin) = prof_1D_node(ibin) + 1
enddo

open(77, file = IO_meshProfile)
do ii = 1, nbin
  write(77,*) ii, prof_1D_node(ii)
enddo
close(77)

! Deallocate memory
deallocate(prof_1D_node)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshProfile
