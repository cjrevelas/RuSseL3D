!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine get_node_volume(volnp, node)
!--------------------------------------------------------------------!
use geometry_mod,     only: ndm, nel, numnp, el_node, global_node_id_type_domain,&
&                       xc, num_of_elems_of_node
use error_handing_mod
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer              :: i, j, ii, l, n, s, lint
integer, intent(out) :: node

real(8), intent(out)         :: volnp
real(8)                      :: xsj, uqp, sumel, vol
real(8), dimension(nel)      :: u_local
real(8), dimension(numnp)    :: u_spat
real(8), dimension(ndm, nel) :: xl
real(8), dimension(4,11)     :: shp
real(8), dimension(5,11)     :: sv
!--------------------------------------------------------------------!
volnp = 0.d0

vol          = 0.d0
u_spat       = 0.d0
u_spat(node) = 1.d0

do s = 1, num_of_elems_of_node(node)
    n = el_node(node, s)

    do i = 1, nel
        ii = global_node_id_type_domain(i,n)

        do j = 1, ndm
            xl(j,i) = xc(j,ii)
        enddo

         u_local(i) = u_spat(ii)
    enddo

    ! Set up for gauss quadrature
    l=3
    call fem_gausspoints(l, lint, sv)

    sumel = 0.d00

    ! Loop over all quadrature points in element
    do l = 1, lint
        call fem_tetshpfun(sv(1,l), xl, ndm, nel, xsj, shp)

        xsj = xsj*sv(5,l)

        uqp = 0.0d0

        do j = 1, nel
             uqp = uqp + shp(4,j)*u_local(j)
        enddo

        sumel = sumel + uqp * xsj
    enddo

    volnp = volnp + sumel
enddo

return
!--------------------------------------------------------------------!
end subroutine get_node_volume
