!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine fem_tetshpfun(xi, xl, ndm, nel, xsj, shp)
!-----------------------------------------------------------------------------------!
use constants_mod, only: tol
!-----------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------!
integer, intent(in) :: ndm, nel
integer             :: nnel, i, j, k

real(8), intent(in), dimension(4)     :: xi
real(8), intent(in), dimension(ndm,*) :: xl
real(8), intent(out), dimension(4,*)  :: shp
real(8), intent(out)                  :: xsj
real(8)                               :: detr
real(8), dimension(3,3)               :: xs, xsi, a
real(8), dimension(3)                 :: te
!-----------------------------------------------------------------------------------!
! Linear shape functions and their derivatives
if(ABS(nel).eq.4) then
     shp(1,1) =  1.d00
     shp(1,2) =  0.d00
     shp(1,3) =  0.d00
     shp(1,4) = -1.d00

     shp(2,1) =  0.d00
     shp(2,2) =  1.d00
     shp(2,3) =  0.d00
     shp(2,4) = -1.d00

     shp(3,1) =  0.d00
     shp(3,2) =  0.d00
     shp(3,3) =  1.d00
     shp(3,4) = -1.d00

     shp(4, 1) = xi(1)
     shp(4, 2) = xi(2)
     shp(4, 3) = xi(3)
     shp(4, 4) = xi(4)

     nnel = 4

! Quadratic shape functions and derivatives
elseif(ABS(nel).eq.10 .or. ABS(nel).eq.11  .or.&
       ABS(nel).eq.14 .or. ABS(nel).eq.15) then

     shp(1, 1) =  4.d0*xi(1) - 1.d0
     shp(1, 2) =  0.d0
     shp(1, 3) =  0.d0
     shp(1, 4) = -4.d0*xi(4) + 1.d0
     shp(1, 5) =  4.d0*xi(2)
     shp(1, 6) =  0.d0
     shp(1, 7) =  4.d0*xi(3)
     shp(1, 8) =  4.d0*(xi(4) - xi(1))
     shp(1, 9) = -4.d0*xi(2)
     shp(1,10) = -4.d0*xi(3)

     shp(2, 1) =  0.d0
     shp(2, 2) =  4.d0*xi(2) - 1.d0
     shp(2, 3) =  0.d0
     shp(2, 4) = -4.d0*xi(4) + 1.d0
     shp(2, 5) =  4.d0*xi(1)
     shp(2, 6) =  4.d0*xi(3)
     shp(2, 7) =  0.d0
     shp(2, 8) = -4.d0*xi(1)
     shp(2, 9) =  4.d0*(xi(4) - xi(2))
     shp(2,10) = -4.d0*xi(3)

     shp(3, 1) =  0.d0
     shp(3, 2) =  0.d0
     shp(3, 3) =  4.d0*xi(3) - 1.d0
     shp(3, 4) = -4.d0*xi(4) + 1.d0
     shp(3, 5) =  0.d0
     shp(3, 6) =  4.d0*xi(2)
     shp(3, 7) =  4.d0*xi(1)
     shp(3, 8) = -4.d0*xi(1)
     shp(3, 9) = -4.d0*xi(2)
     shp(3,10) =  4.d0*(xi(4) - xi(3))

     ! Compute shape functions and store in shp(4,*)
     shp(4, 1) = (2.d0*xi(1) - 1.d0) * xi(1)
     shp(4, 2) = (2.d0*xi(2) - 1.d0) * xi(2)
     shp(4, 3) = (2.d0*xi(3) - 1.d0) * xi(3)
     shp(4, 4) = (2.d0*xi(4) - 1.d0) * xi(4)
     shp(4, 5) = 4.d0 * xi(1) * xi(2)
     shp(4, 6) = 4.d0 * xi(2) * xi(3)
     shp(4, 7) = 4.d0 * xi(3) * xi(1)
     shp(4, 8) = 4.d0 * xi(4) * xi(1)
     shp(4, 9) = 4.d0 * xi(4) * xi(2)
     shp(4,10) = 4.d0 * xi(4) * xi(3)

     nnel = nel
else
    write(*,2000) nel
endif

! Compute jacobian matrix
do i = 1,3
   do j = 1,3
        xs(j,i) = 0.0d0
        do k = 1,nnel
            xs(j,i) = xs(j,i) - xl(j,k)*shp(i,k)
        enddo
   enddo
enddo

! Compute inverse of jacobian matrix
xsi(1,1) = xs(2,2)*xs(3,3) - xs(3,2)*xs(2,3)
xsi(1,2) = xs(3,2)*xs(1,3) - xs(1,2)*xs(3,3)
xsi(1,3) = xs(1,2)*xs(2,3) - xs(2,2)*xs(1,3)

xsi(2,1) = xs(2,3)*xs(3,1) - xs(3,3)*xs(2,1)
xsi(2,2) = xs(3,3)*xs(1,1) - xs(1,3)*xs(3,1)
xsi(2,3) = xs(1,3)*xs(2,1) - xs(2,3)*xs(1,1)

xsi(3,1) = xs(2,1)*xs(3,2) - xs(3,1)*xs(2,2)
xsi(3,2) = xs(3,1)*xs(1,2) - xs(1,1)*xs(3,2)
xsi(3,3) = xs(1,1)*xs(2,2) - xs(2,1)*xs(1,2)

do i = 1,3
   do j = 1,3
      a(i,j) = xs(i,1)*xsi(1,j)+xs(i,2)*xsi(2,j)+xs(i,3)*xsi(3,j)
   enddo
enddo

! Compute jacobian determinant
xsj = xs(1,1)*xsi(1,1) + xs(1,2)*xsi(2,1) + xs(1,3)*xsi(3,1)

if (DABS(xsj)>tol) then
    detr = 1.d0/xsj
    xsj  = xsj!*(1./6.)
else
    !write(iow,*) ' TETSHP: Determinant =',xsj
    detr = 1.d0
endif

! Compute jacobian inverse
do j = 1, 3
   do i = 1, 3
        xs(i,j) = xsi(i,j)*detr
   enddo
enddo

! Compute shape function derivatives wrt to real node coordinates
do k = 1, nnel
   do i = 1, 3
        te(i) = shp(1,k)*xs(1,i) + shp(2,k)*xs(2,i) + shp(3,k)*xs(3,i)
   enddo

   do i = 1,3
        shp(i,k) = te(i)
   enddo
enddo

2000  format(/" *ERROR* TETSHP not coded for nel =",I4)
!-----------------------------------------------------------------------------------!
end subroutine fem_tetshpfun
