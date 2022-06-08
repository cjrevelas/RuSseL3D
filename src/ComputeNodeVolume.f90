!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ComputeNodeVolume(volnp, node)
!--------------------------------------------------------------------!
use geometry_mod,     only: numDimensions, numNodesLocalTypeDomain, &
                            numNodes, elementOfNode, nodeCoord,     &
                            globalNodeIdTypeDomain, numElementsOfNode
use error_handing_mod
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer              :: ii, jj, kk, ll, nn, ss, lint
integer, intent(out) :: node

real(8), intent(out)                                       :: volnp
real(8), dimension(numNodesLocalTypeDomain)                :: u_local
real(8), dimension(numNodes)                               :: u_spat
real(8), dimension(numDimensions, numNodesLocalTypeDomain) :: xl
real(8), dimension(4,11)                                   :: shp
real(8), dimension(5,11)                                   :: sv
real(8)                                                    :: xsj, uqp
real(8)                                                    :: sumel, vol
!--------------------------------------------------------------------!
volnp = 0.0d0

vol          = 0.0d0
u_spat       = 0.0d0
u_spat(node) = 1.0d0

do ss = 1, numElementsOfNode(node)
  nn = elementOfNode(node, ss)

  do kk = 1, numNodesLocalTypeDomain
    ii = globalNodeIdTypeDomain(kk,nn)

    do jj = 1, numDimensions
      xl(jj,kk) = nodeCoord(jj,ii)
    enddo

    u_local(kk) = u_spat(ii)
  enddo

  ! Set up for gauss quadrature
  ll = 3
  call FemGaussPoints(ll, lint, sv)

  sumel = 0.0d0

  ! Loop over all quadrature points in element
  do ll = 1, lint
    call FemShapeFunctions(sv(1,ll), xl, numDimensions, numNodesLocalTypeDomain, xsj, shp)

    xsj = xsj*sv(5,ll)

    uqp = 0.0d0

    do jj = 1, numNodesLocalTypeDomain
      uqp = uqp + shp(4,jj)*u_local(jj)
    enddo

    sumel = sumel + uqp * xsj
  enddo

  volnp = volnp + sumel
enddo

return
!--------------------------------------------------------------------!
end subroutine ComputeNodeVolume
