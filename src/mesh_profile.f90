!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine mesh_profile()
!----------------------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod, only: prof_dim
use geometry_mod,    only: numnp, box_lo, box_hi, xc
use iofiles_mod,     only: mesh_prof
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
integer :: ii, ibin, nbin
integer, allocatable, dimension(:) :: prof_1D_node

real(8) :: prof_bin
!----------------------------------------------------------------------------------------------------------------------------------!
prof_bin = 0.5d0
nbin     = NINT((box_hi(prof_dim) - box_lo(prof_dim)) / prof_bin) + 1

allocate(prof_1D_node(nbin))

prof_1D_node=0
do ii = 1, numnp
    ibin               = NINT((xc(prof_dim,ii) - box_lo(prof_dim))/prof_bin) + 1
    prof_1D_node(ibin) = prof_1D_node(ibin) + 1
enddo

open(77, file = mesh_prof)
do ii = 1, nbin
    write(77,*) ii, prof_1D_node(ii)
enddo
close(77)

! Deallocate memory
deallocate(prof_1D_node)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine mesh_profile
