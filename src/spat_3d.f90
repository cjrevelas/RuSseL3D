subroutine spat_3d(u_spat, sum_out, Q)
!--------------------------------------------------------------------!
use xdata
use kcw
use constants
use error_handing
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: lint
integer :: i, j, ii, l, n

real(8), intent(in), dimension(numnp) :: u_spat
real(8), intent(out)                  :: sum_out, Q
real(8)                               :: xsj, uqp, sumel, volel, vol
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
    call gauss_3d(l, lint, sv)

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

Q = sum_out/vol

!write(iow,'('' SUM '',E16.9,'' volume '',E16.9,'' mean value  of solution '',E16.9)') sum_out, vol, Q
!write(*,'('' SUM '',E16.9,'' volume '',E16.9,'' mean value  of solution '',E16.9)') sum_out, vol, Q

if (dabs(vol - volume) > 1e-6) then
    write(ERROR_MESSAGE,'(''Warning! Spat volume ('',E16.9,'') is different than mesh volume ('',E16.9,'')'')')vol, volume
    call exit_with_error(0,2,0,ERROR_MESSAGE)
endif


return
!--------------------------------------------------------------------!
end subroutine spat_3d
