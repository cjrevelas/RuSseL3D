subroutine get_node_volume(volnp, node)
!--------------------------------------------------------------------!
use geometry,     only: ndm, nel, numnp, el_node, ix, xc, n_el_node
use error_handing
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer              :: i, j, ii, l, n, s, lint
integer, intent(out) :: node

real(8), intent(out)         :: volnp
real(8)                      :: xsj, uqp, sumel, vol
real(8), dimension(nel)      :: ul
real(8), dimension(numnp)    :: u_spat
real(8), dimension(ndm, nel) :: xl
real(8), dimension(4,11)     :: shp
real(8), dimension(5,11)     :: sv
!--------------------------------------------------------------------!
!initialize accumulator for integral of the solution
volnp = 0.d0

!initialize accumulator for total domain volume
vol          = 0.d0
u_spat       = 0.d0
u_spat(node) = 1.d0

do s = 1, n_el_node(node)
    n = el_node(node, s)

    do i = 1, nel
        !find global index of node
        ii = ix(i,n)

        do j = 1, ndm
            xl(j,i) = xc(j,ii)
        enddo

        !copy value of solution at current node from global to local array
         ul(i) = u_spat(ii)
    enddo

    !set up for gauss quadrature
    l=3
    call fem_gausspoints(l, lint, sv)

    sumel = 0.d00

    !loop over all quadrature points in element
    do l = 1, lint
        call fem_tetshpfun(sv(1,l), xl, ndm, nel, xsj, shp)

        xsj = xsj*sv(5,l)

        uqp = 0.0d0

        do j = 1, nel
             uqp = uqp + shp(4,j)*ul(j)
        enddo

        sumel = sumel + uqp * xsj
    enddo

    volnp = volnp + sumel
enddo

return
!--------------------------------------------------------------------!
end subroutine get_node_volume
