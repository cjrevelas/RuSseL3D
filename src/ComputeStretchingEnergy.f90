!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

real(8) function ComputeStretchingEnergy(gnode_id, qmx, qgr)
!-------------------------------------------------------------------------------------------------!
use constants_mod,   only: A3_to_m3
use parser_vars_mod, only: beta, numConvolPointsGrafted, lengthGrafted, &
                           rg2OfGraftedMonomer, segmentBulkDensity
use geometry_mod,    only: numNodes, nodeCoord
!-------------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------------!
integer, intent(in) :: gnode_id
integer             :: ii

real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qmx, qgr
real(8), dimension(numNodes)                                      :: dx2, dy2, dz2, dr2
real(8), dimension(numNodes)                                      :: phi_end, rho_end, A_stretch
real(8)                                                           :: Q, vol
!-------------------------------------------------------------------------------------------------!
do ii = 1, numNodes
  phi_end(ii) = 1 / lengthGrafted * qgr(numConvolPointsGrafted+1,ii) * qmx(1,ii)
  rho_end(ii) = phi_end(ii) * segmentBulkDensity * A3_to_m3
  dx2(ii)     = (nodeCoord(1,ii)-nodeCoord(1,gnode_id))**2.0d0
  dy2(ii)     = (nodeCoord(2,ii)-nodeCoord(2,gnode_id))**2.0d0
  dz2(ii)     = (nodeCoord(3,ii)-nodeCoord(3,gnode_id))**2.0d0
  dr2(ii)     = dx2(ii) + dy2(ii) + dz2(ii)

  A_stretch(ii) = 3.0d0/(2.0d0 * beta * (rg2OfGraftedMonomer*lengthGrafted*6.0d0)) * dr2(ii)
enddo

call FemIntegration(rho_end*A_stretch, ComputeStretchingEnergy, Q, vol)
!-------------------------------------------------------------------------------------------------!
end function ComputeStretchingEnergy
