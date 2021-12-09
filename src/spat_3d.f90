subroutine spat_3d(u_spat, sum_out, QQ, vol)
!--------------------------------------------------------------------!
use geometry, only : numnp, ndm, nel, numel, ix, xc
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: lint
integer :: i, j, ii, l, n

real(8), intent(in), dimension(numnp) :: u_spat
real(8), intent(out)                  :: sum_out, QQ, vol
real(8)                               :: xsj, uqp, sumel, volel
real(8), dimension(nel)               :: ul
real(8), dimension(ndm, nel)          :: xl
real(8), dimension(4,11)              :: shp
real(8), dimension(5,11)              :: sv
!--------------------------------------------------------------------!
!initialize accumulator for integral of the solution
sum_out = 0.d00

!initialize accumulator for total domain volume
vol = 0.d00

do n = 1, numel

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
    call gausspoints(l, lint, sv)

    sumel = 0.d00
    volel = 0.d00

    !loop over all quadrature points in element
    do l = 1, lint

        call tetshp(sv(1,l), xl, ndm, nel, xsj, shp)

        xsj = xsj*sv(5,l)

        uqp = 0.0d0

        do j = 1, nel
             uqp = uqp + shp(4,j)*ul(j)
        enddo

        sumel = sumel + uqp * xsj

        volel = volel + xsj

    enddo

    sum_out = sum_out + sumel

    vol = vol + volel

enddo !n

QQ = sum_out/vol

return
!--------------------------------------------------------------------!
end subroutine spat_3d
