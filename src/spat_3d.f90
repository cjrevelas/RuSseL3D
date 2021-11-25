subroutine spat_3d(u_spat, sum_out, Q)
!--------------------------------------------------------------------!
use xdata
use mdata
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: lint

integer :: i, j, ii, l, n

real(8) :: Q, xsj

real(8), dimension(nel) :: ul
real(8), dimension(ndm, nel) :: xl

real(8), dimension(4,11) :: shp
real(8), dimension(5,11) :: sv
real(8), dimension(numnp) :: u_spat

real(8) :: uqp, sumel, sum_out, volel, vol
!--------------------------------------------------------------------!

!Initialize accumulator for integral of the solution
sum_out = 0.d00

!Initialize accumulator for total domain volume
vol = 0.d00


do n = 1, numel

    do i = 1, nel
        !Find global index of node
        ii = ix(i,n)

        do j = 1, ndm
            xl(j,i) = xc(j,ii)
        enddo

        !Copy value of solution at current node from global to local array
         ul(i) = u_spat(ii)
    enddo

    !Set up for gauss quadrature.
    l=3
    call gauss_3d(l, lint, sv)

    sumel = 0.d00
    volel = 0.d00

    !Loop over all quadrature points in element
    do l = 1, lint

        !call tetshp10(sv(1,l), xl, ndm, xsj, shp)
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

Q = sum_out/vol

write(iow,'('' SUM '',E16.9,'' volume '',E16.9,'' mean value  of solution '',E16.9)') sum_out, vol, q

return
!--------------------------------------------------------------------!
end subroutine spat_3d
