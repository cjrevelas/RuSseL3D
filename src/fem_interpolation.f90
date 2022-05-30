!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

function fem_interpolation(nodeId, x_interp, y_interp, z_interp, uu) result(u_interp)
!--------------------------------------------------------------------------------------------------------------------------------------------!
use geometry_mod, only : numNodes, numDimensions, numNodesLocalTypeDomain, elementOfNode, globalNodeIdTypeDomain, nodeCoord, numElementsOfNode
!--------------------------------------------------------------------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------------------------------------------------------------------!
integer, intent(in)                   :: nodeId
integer                               :: jj, ll, lint, nn, m1, m2, inside, mm, iglob, iloc

real(8), intent(in), dimension(numNodes) :: uu
real(8), intent(in)                      :: x_interp, y_interp, z_interp
real(8), dimension(4,11)                 :: shp
real(8), dimension(5,11)                 :: sv
real(8), dimension(3,4)                  :: xl
real(8), dimension(4)                    :: ul_uu
real(8), dimension(4)                    :: gc, lc
real(8), dimension(4,4)                  :: transf, transf_inv
real(8)                                  :: volel, xsj, u_interp
!--------------------------------------------------------------------------------------------------------------------------------------------!
interface
    pure function tools_matinv4(transf)
    real(8), intent(in), dimension(4,4) :: transf
    real(8),             dimension(4,4) :: tools_matinv4
    end function
end interface
!--------------------------------------------------------------------------------------------------------------------------------------------!
ll = 3
call fem_gausspoints(ll, lint, sv)

u_interp = 0.0d0

do mm = 1, numElementsOfNode(nodeId)
  nn = elementOfNode(nodeId, mm)

  do iloc = 1, 4
    iglob = globalNodeIdTypeDomain(iloc,nn)

    do jj = 1, 3
      xl(jj,iloc) = nodeCoord(jj,iglob)
    enddo

    ul_uu(iloc)  = uu(iglob)
  enddo

  volel = 0.0d0

  do ll = 1, lint
    call fem_tetshpfun(sv(1,ll), xl, numDimensions, numNodesLocalTypeDomain, xsj, shp)

    xsj = xsj*sv(5,ll)

    volel = volel + xsj
  enddo

  call fem_is_node_inside_el(xl, x_interp, y_interp, z_interp, lint, sv, volel, inside, shp, numDimensions, numNodesLocalTypeDomain)

  if (inside==1) then
    gc(1) = 1.0d0
    gc(2) = x_interp
    gc(3) = y_interp
    gc(4) = z_interp

    do m1 = 1, 4
      transf(1,m1) = 1.0d0
      transf(2,m1) = xl(1,m1)
      transf(3,m1) = xl(2,m1)
      transf(4,m1) = xl(3,m1)
    enddo

    transf_inv = tools_matinv4(transf)

    do m1 = 1, 4
      lc(m1) = 0.0d0
      do m2 = 1, 4
        lc(m1) = lc(m1) + transf_inv(m1,m2) * gc(m2)
      enddo
    enddo

    call fem_tetshpfun(lc, xl, numDimensions, numNodesLocalTypeDomain, xsj, shp)

    u_interp  = 0.0d0

    do jj = 1, 4
      u_interp  = u_interp  + shp(4,jj) * ul_uu(jj)
    enddo

    exit
  endif
enddo
!--------------------------------------------------------------------------------------------------------------------------------------------!
end function fem_interpolation
