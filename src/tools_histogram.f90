!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine tools_histogram(upd_lbin, volnp)
!-----------------------------------------------------------------------------------------------------------!
use hist_mod
use geometry_mod,    only: numNodes, boxLength, isDirichletFace, boxLow, boxHigh, nodeCoord
use parser_vars_mod, only: numNanoparticleFaces, wall_distance, center_np, radius_np_eff
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: kk, mm, nn, bin

real(8), intent(in)                   :: upd_lbin
real(8), intent(in), dimension(numNodes) :: volnp
real(8)                               :: r_center_surf, r_centers, radius_np_actual, Rmax
!-----------------------------------------------------------------------------------------------------------!
lbin = upd_lbin
Rmax = SQRT(boxLength(1)**2.0d0 + boxLength(2)**2.0d0 + boxLength(3)**2.0d0)
nbin = INT(Rmax / lbin) + 1

if (ALLOCATED(planar_cell_of_np)) deallocate(planar_cell_of_np)
if (ALLOCATED(sph_cell_of_np))    deallocate(sph_cell_of_np)
if (ALLOCATED(dist_from_face))    deallocate(dist_from_face)
if (ALLOCATED(dist_from_np))      deallocate(dist_from_np)
if (ALLOCATED(cell_vol_planar))   deallocate(cell_vol_planar)
if (ALLOCATED(cell_vol_sph))      deallocate(cell_vol_sph)

allocate(planar_cell_of_np(numNodes,3,2))
allocate(dist_from_face(numNodes,3,2))
allocate(cell_vol_planar(nbin,3,2))
allocate(sph_cell_of_np(numNanoparticleFaces,numNodes))
allocate(dist_from_np(numNanoparticleFaces,numNodes))
allocate(cell_vol_sph(numNanoparticleFaces,nbin))

planar_cell_of_np = 0
sph_cell_of_np    = 0
dist_from_face    = 0.0d0
dist_from_np      = 0.0d0
cell_vol_planar   = 0.0d0
cell_vol_sph      = 0.0d0
r_center_surf     = 0.0d0

! Binning in planar geometries
do mm = 1, 3
  do nn = 1, 2
    if (isDirichletFace(mm,nn)) then
      do kk = 1, numNodes
        if (nn.eq.1) then
          r_center_surf = nodeCoord(mm,kk) - boxLow(mm) + wall_distance
        elseif (nn.eq.2) then
          r_center_surf = boxHigh(mm) - nodeCoord(mm,kk) + wall_distance
        endif
        bin                         = INT(r_center_surf/lbin)+1
        planar_cell_of_np(kk,mm,nn) = bin
        cell_vol_planar(bin,mm,nn)  = cell_vol_planar(bin,mm,nn) + volnp(kk)
        dist_from_face(kk,mm,nn)    = r_center_surf
      enddo
    endif
  enddo
enddo

! Binning in spherical geometries
do mm = 1, numNanoparticleFaces
  do kk = 1, numNodes
    r_centers             = DSQRT((nodeCoord(1,kk)-center_np(1,mm))**2.0d0 + (nodeCoord(2,kk)-center_np(2,mm))**2.0d0 &
                                                                           + (nodeCoord(3,kk)-center_np(3,mm))**2.0d0)
    radius_np_actual      = radius_np_eff(mm) - wall_distance
    r_center_surf         = r_centers - radius_np_actual
    bin                   = INT(r_center_surf/lbin)+1
    sph_cell_of_np(mm,kk) = bin
    cell_vol_sph(mm,bin)  = cell_vol_sph(mm,bin) + volnp(kk)
    dist_from_np(mm,kk)   = r_center_surf
  enddo
enddo
!-----------------------------------------------------------------------------------------------------------!
end subroutine tools_histogram
