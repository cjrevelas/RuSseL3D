subroutine tetshp(xi, xl, xsj, shp)
!--------------------------------------------------------------------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------------------------------------------------------------------! 
integer :: i, j, k

real(8), intent(inout), dimension(4) :: xi
real(8), intent(in), dimension(3,*)  :: xl
real(8), intent(out), dimension(4,*) :: shp
real(8), intent(out)                 :: xsj
real(8), dimension(3,3)              :: xs, xsinv
real(8), dimension(3,10)             :: shpxi
!--------------------------------------------------------------------------------------------------------------------------------------------! 
!in case that we pass in something different than the gauss points
xi(4) = 1.d0 - xi(1) - xi(2) - xi(3)

shp(4,1) = xi(1)
shp(4,2) = xi(2)
shp(4,3) = xi(3)
shp(4,4) = xi(4)

shpxi(1,1) =  1.d00
shpxi(1,2) =  0.d00
shpxi(1,3) =  0.d00
shpxi(1,4) = -1.d00

shpxi(2,1) =  0.d00
shpxi(2,2) =  1.d00
shpxi(2,3) =  0.d00
shpxi(2,4) = -1.d00

shpxi(3,1) =  0.d00
shpxi(3,2) =  0.d00
shpxi(3,3) =  1.d00
shpxi(3,4) = -1.d00

!compute transpose of Jacobian matrix of the isoparametric transformation
!
!                  partial r_j    10          partial N^k
!xs(i,j) = J_ji = ------------- = Sum xl(j,k) ------------
!                  partial xi_i   k=1         partial xi_i
!
!i,j in {1,2,3}

do i = 1, 3
    do j = 1, 3
        xs(i,j) = 0.0d0
        do k = 1, 4
            xs(i,j) = xs(i,j) - shpxi(i,k)*xl(j,k)
        enddo
    enddo
enddo

!compute jacobian determinant
xsj  = xs(1,1)*(xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2)) &
     + xs(1,2)*(xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)) &
     + xs(1,3)*(xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1))

!compute Jacobian inverse:
!xsinv(i,j) = (J^-1)_ij = partial xi_i/partial r_j

xsinv(1,1) = xs(2,2)*xs(3,3) - xs(3,2)*xs(2,3)
xsinv(2,1) = xs(3,2)*xs(1,3) - xs(1,2)*xs(3,3)
xsinv(3,1) = xs(1,2)*xs(2,3) - xs(2,2)*xs(1,3)

xsinv(1,2) = xs(3,1)*xs(2,3) - xs(2,1)*xs(3,3)
xsinv(2,2) = xs(1,1)*xs(3,3) - xs(3,1)*xs(1,3)
xsinv(3,2) = xs(2,1)*xs(1,3) - xs(1,1)*xs(2,3)

xsinv(1,3) = xs(2,1)*xs(3,2) - xs(3,1)*xs(2,2)
xsinv(2,3) = xs(3,1)*xs(1,2) - xs(1,1)*xs(3,2)
xsinv(3,3) = xs(1,1)*xs(2,2) - xs(2,1)*xs(1,2)

do j = 1, 3
   do i = 1, 3
      xsinv(i,j) = xsinv(i,j)/xsj
   enddo
enddo

!compute derivatives of shape functions with respect to nodal
!coordinates.

!shp(j,k) = partial N^k/partial r_j     (j=1,2,3)
!
!here it is calculated as

!              3   partial N^k   partial xi_i     3
!  shp(j,k) = Sum ------------- -------------- = Sum shpxi(i,k) xsinv(i,j)
!             i=1  partial xi_i  partial r_j     i=1

do k = 1, 4
    do j = 1, 3
        shp(j,k) = 0.0d0
        do i = 1,3
            shp(j,k) = shp(j,k)+shpxi(i,k)*xsinv(i,j)
        enddo
    enddo
enddo
!--------------------------------------------------------------------------------------------------------------------------------------------!
end subroutine tetshp
