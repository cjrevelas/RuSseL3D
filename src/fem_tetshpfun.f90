!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine fem_tetshpfun(xi, xl, numDimensions, nel, xsj, shp)
!-----------------------------------------------------------------------------------!
use constants_mod, only: tol
!-----------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------!
integer, intent(in) :: numDimensions, nel
integer             :: nnel, ii, jj, kk

real(8), intent(in), dimension(4)               :: xi
real(8), intent(in), dimension(numDimensions,*) :: xl
real(8), intent(out), dimension(4,*)            :: shp
real(8), intent(out)                            :: xsj
real(8)                                         :: detr
real(8), dimension(3,3)                         :: xs, xsi, a
real(8), dimension(3)                           :: te
!-----------------------------------------------------------------------------------!
! Linear shape functions and their derivatives
if(ABS(nel).eq.4) then
  shp(1,1) =  1.0d0
  shp(1,2) =  0.0d0
  shp(1,3) =  0.0d0
  shp(1,4) = -1.0d0

  shp(2,1) =  0.0d0
  shp(2,2) =  1.0d0
  shp(2,3) =  0.0d0
  shp(2,4) = -1.0d0

  shp(3,1) =  0.0d0
  shp(3,2) =  0.0d0
  shp(3,3) =  1.0d0
  shp(3,4) = -1.0d0

  shp(4, 1) = xi(1)
  shp(4, 2) = xi(2)
  shp(4, 3) = xi(3)
  shp(4, 4) = xi(4)

  nnel = 4

! Quadratic shape functions and derivatives
elseif(ABS(nel).eq.10 .or. ABS(nel).eq.11 .or. ABS(nel).eq.14 .or. ABS(nel).eq.15) then
  shp(1, 1) =  4.0d0*xi(1) - 1.0d0
  shp(1, 2) =  0.0d0
  shp(1, 3) =  0.0d0
  shp(1, 4) = -4.0d0*xi(4) + 1.0d0
  shp(1, 5) =  4.0d0*xi(2)
  shp(1, 6) =  0.0d0
  shp(1, 7) =  4.0d0*xi(3)
  shp(1, 8) =  4.0d0*(xi(4) - xi(1))
  shp(1, 9) = -4.0d0*xi(2)
  shp(1,10) = -4.0d0*xi(3)

  shp(2, 1) =  0.0d0
  shp(2, 2) =  4.0d0*xi(2) - 1.0d0
  shp(2, 3) =  0.0d0
  shp(2, 4) = -4.0d0*xi(4) + 1.0d0
  shp(2, 5) =  4.0d0*xi(1)
  shp(2, 6) =  4.0d0*xi(3)
  shp(2, 7) =  0.0d0
  shp(2, 8) = -4.0d0*xi(1)
  shp(2, 9) =  4.0d0*(xi(4) - xi(2))
  shp(2,10) = -4.0d0*xi(3)

  shp(3, 1) =  0.0d0
  shp(3, 2) =  0.0d0
  shp(3, 3) =  4.0d0*xi(3) - 1.0d0
  shp(3, 4) = -4.0d0*xi(4) + 1.0d0
  shp(3, 5) =  0.0d0
  shp(3, 6) =  4.0d0*xi(2)
  shp(3, 7) =  4.0d0*xi(1)
  shp(3, 8) = -4.0d0*xi(1)
  shp(3, 9) = -4.0d0*xi(2)
  shp(3,10) =  4.0d0*(xi(4) - xi(3))

  ! Compute shape functions and store in shp(4,*)
  shp(4, 1) = (2.0d0*xi(1) - 1.0d0) * xi(1)
  shp(4, 2) = (2.0d0*xi(2) - 1.0d0) * xi(2)
  shp(4, 3) = (2.0d0*xi(3) - 1.0d0) * xi(3)
  shp(4, 4) = (2.0d0*xi(4) - 1.0d0) * xi(4)
  shp(4, 5) = 4.0d0 * xi(1) * xi(2)
  shp(4, 6) = 4.0d0 * xi(2) * xi(3)
  shp(4, 7) = 4.0d0 * xi(3) * xi(1)
  shp(4, 8) = 4.0d0 * xi(4) * xi(1)
  shp(4, 9) = 4.0d0 * xi(4) * xi(2)
  shp(4,10) = 4.0d0 * xi(4) * xi(3)

  nnel = nel
else
  write(*,2000) nel
endif

! Compute jacobian matrix
do ii = 1, 3
  do jj = 1, 3
    xs(jj,ii) = 0.0d0
    do kk = 1, nnel
      xs(jj,ii) = xs(jj,ii) - xl(jj,kk)*shp(ii,kk)
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

do ii = 1, 3
  do jj = 1, 3
    a(ii,jj) = xs(ii,1)*xsi(1,jj)+xs(ii,2)*xsi(2,jj)+xs(ii,3)*xsi(3,jj)
  enddo
enddo

! Compute jacobian determinant
xsj = xs(1,1)*xsi(1,1) + xs(1,2)*xsi(2,1) + xs(1,3)*xsi(3,1)

if (DABS(xsj)>tol) then
  detr = 1.0d0/xsj
  xsj  = xsj!*(1.0d0/6.0d0)
else
  !write(iow,*) ' TETSHP: Determinant =',xsj
  detr = 1.0d0
endif

! Compute jacobian inverse
do jj = 1, 3
  do ii = 1, 3
    xs(ii,jj) = xsi(ii,jj)*detr
  enddo
enddo

! Compute shape function derivatives wrt to real node coordinates
do kk = 1, nnel
  do ii = 1, 3
    te(ii) = shp(1,kk)*xs(1,ii) + shp(2,kk)*xs(2,ii) + shp(3,kk)*xs(3,ii)
  enddo

  do ii = 1,3
    shp(ii,kk) = te(ii)
  enddo
enddo

2000  format(/" *ERROR* TETSHP not coded for nel =",I4)
!-----------------------------------------------------------------------------------!
end subroutine fem_tetshpfun
