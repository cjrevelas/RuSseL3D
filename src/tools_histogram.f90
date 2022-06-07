!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine tools_histogram(upd_lbin, volnp)
!-----------------------------------------------------------------------------------------------------------!
use hist_mod
use geometry_mod,    only: numNodes, boxLength, isDirichletFace, boxLow, boxHigh, nodeCoord
use parser_vars_mod, only: numNanopFaces, wallDistance, nanopCenter, nanopRadiusEff
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: kk, mm, nn, bin

real(8), intent(in), dimension(numNodes) :: volnp
real(8), intent(in)                      :: upd_lbin
real(8)                                  :: r_center_surf, r_centers, nanopRadiusActual, Rmax
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
allocate(sph_cell_of_np(numNanopFaces,numNodes))
allocate(dist_from_np(numNanopFaces,numNodes))
allocate(cell_vol_sph(numNanopFaces,nbin))

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
          r_center_surf = nodeCoord(mm,kk) - boxLow(mm) + wallDistance
        elseif (nn.eq.2) then
          r_center_surf = boxHigh(mm) - nodeCoord(mm,kk) + wallDistance
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
do mm = 1, numNanopFaces
  do kk = 1, numNodes
    r_centers             = DSQRT((nodeCoord(1,kk)-nanopCenter(1,mm))**2.0d0 + (nodeCoord(2,kk)-nanopCenter(2,mm))**2.0d0 &
                                                                             + (nodeCoord(3,kk)-nanopCenter(3,mm))**2.0d0)
    nanopRadiusActual     = nanopRadiusEff(mm) - wallDistance
    r_center_surf         = r_centers - nanopRadiusActual
    bin                   = INT(r_center_surf/lbin)+1
    sph_cell_of_np(mm,kk) = bin
    cell_vol_sph(mm,bin)  = cell_vol_sph(mm,bin) + volnp(kk)
    dist_from_np(mm,kk)   = r_center_surf
  enddo
enddo
!-----------------------------------------------------------------------------------------------------------!
end subroutine tools_histogram
