!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine fem_integration(u_spat, sum_out, QQ, vol)
!--------------------------------------------------------------------!
use geometry_mod, only : numnp, ndm, nel, numel, global_node_id_type_domain, xc
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: lint
integer :: i, j, ii, l, n

real(8), intent(in), dimension(numnp) :: u_spat
real(8), intent(out)                  :: sum_out, QQ, vol
real(8)                               :: xsj, uqp, sumel, volel
real(8), dimension(nel)               :: u_local
real(8), dimension(ndm, nel)          :: xl
real(8), dimension(4,11)              :: shp
real(8), dimension(5,11)              :: sv
!--------------------------------------------------------------------!
sum_out = 0.0d0
vol     = 0.0d0

do n = 1, numel

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
    volel = 0.d00

    ! Loop over all quadrature points in element
    do l = 1, lint

        call fem_tetshpfun(sv(1,l), xl, ndm, nel, xsj, shp)

        xsj = xsj*sv(5,l)

        uqp = 0.0d0

        do j = 1, nel
             uqp = uqp + shp(4,j)*u_local(j)
        enddo

        sumel = sumel + uqp * xsj

        volel = volel + xsj

    enddo

    sum_out = sum_out + sumel

    vol = vol + volel

enddo

QQ = sum_out/vol

return
!--------------------------------------------------------------------!
end subroutine fem_integration
