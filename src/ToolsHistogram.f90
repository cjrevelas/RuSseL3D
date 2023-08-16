!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ToolsHistogram(updatedBinLength, nodeVolume)
!-----------------------------------------------------------------------------------------------------------!
use hist_mod
use geometry_mod,    only: numNodes, boxLength, isDirichletFace, boxLow, boxHigh, nodeCoord
use parser_vars_mod, only: numNanopFaces, wallDistance, nanopCenter, nanopRadiusEff
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: axis, face, node, bin

real(8), intent(in), dimension(numNodes) :: nodeVolume
real(8), intent(in)                      :: updatedBinLength
real(8)                                  :: nanopRadiusActual, maxDistance
real(8)                                  :: centerDistance, surfaceDistance
!-----------------------------------------------------------------------------------------------------------!
binLength   = updatedBinLength
maxDistance = SQRT(boxLength(1)**2.0d0 + boxLength(2)**2.0d0 + boxLength(3)**2.0d0)
numBins     = INT(maxDistance / binLength) + 1

if (ALLOCATED(planarCellId))        deallocate(planarCellId)
if (ALLOCATED(sphericalCellId))     deallocate(sphericalCellId)
if (ALLOCATED(distanceFromFace))    deallocate(distanceFromFace)
if (ALLOCATED(distanceFromNanop))   deallocate(distanceFromNanop)
if (ALLOCATED(planarCellVolume))    deallocate(planarCellVolume)
if (ALLOCATED(sphericalCellVolume)) deallocate(sphericalCellVolume)

allocate(planarCellId(numNodes,3,2))
allocate(distanceFromFace(numNodes,3,2))
allocate(planarCellVolume(numBins,3,2))
allocate(sphericalCellId(numNanopFaces,numNodes))
allocate(distanceFromNanop(numNanopFaces,numNodes))
allocate(sphericalCellVolume(numNanopFaces,numBins))

planarCellId        = 0
sphericalCellId     = 0
distanceFromFace    = 0.0d0
distanceFromNanop   = 0.0d0
planarCellVolume    = 0.0d0
sphericalCellVolume = 0.0d0
surfaceDistance     = 0.0d0

! Binning in planar geometries
do axis = 1, 3
  do face = 1, 2
    if (isDirichletFace(axis,face)) then
      do node = 1, numNodes
        if (face.eq.1) then
          surfaceDistance = nodeCoord(axis,node) - boxLow(axis) + wallDistance
        elseif (face.eq.2) then
          surfaceDistance = boxHigh(axis) - nodeCoord(axis,node) + wallDistance
        endif
        bin                              = INT(surfaceDistance/binLength) + 1
        planarCellId(node,axis,face)     = bin
        planarCellVolume(bin,axis,face)  = planarCellVolume(bin,axis,face) + nodeVolume(node)
        distanceFromFace(node,axis,face) = surfaceDistance
      enddo
    endif
  enddo
enddo

! Binning in spherical geometries
do face = 1, numNanopFaces
  do node = 1, numNodes
    centerDistance                = DSQRT((nodeCoord(1,node)-nanopCenter(1,face))**2.0d0 + (nodeCoord(2,node)-nanopCenter(2,face))**2.0d0 &
                                                                                         + (nodeCoord(3,node)-nanopCenter(3,face))**2.0d0)
    nanopRadiusActual             = nanopRadiusEff(face) - wallDistance
    surfaceDistance               = centerDistance - nanopRadiusActual
    bin                           = INT(surfaceDistance/binLength) + 1
    sphericalCellId(face,node)    = bin
    sphericalCellVolume(face,bin) = sphericalCellVolume(face,bin) + nodeVolume(node)
    distanceFromNanop(face,node)  = surfaceDistance
  enddo
enddo
!-----------------------------------------------------------------------------------------------------------!
end subroutine ToolsHistogram
