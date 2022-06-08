!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine FemIsPointInElement(xl, x_test, y_test, z_test, lint, sv, volel, inside, shp, numDimensions, nel)
!--------------------------------------------------------------------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------------------------------------------------------------------!
integer, intent(in)  :: lint, numDimensions, nel
integer, intent(out) :: inside
integer              :: pp, kk, ll

real(8), intent(in)                    :: x_test, y_test, z_test
real(8), intent(in), dimension(5,11)   :: sv
real(8), intent(inout), dimension(3,4) :: xl
real(8), intent(out), dimension(4,11)  :: shp
real(8), intent(out)                   :: volel
real(8)                                :: temp_x, temp_y, temp_z, vol_test, xsj
!--------------------------------------------------------------------------------------------------------------------------------------------!
! Replace each coordinate vertex with the test point and calculate the new volume each time
inside = 0
kk     = 0

!open(unit = 60, file = "d.elem_vol")
do pp = 1, 4
  vol_test = 0.0d0

  temp_x   = xl(1,pp)
  xl(1,pp) = x_test

  temp_y   = xl(2,pp)
  xl(2,pp) = y_test

  temp_z   = xl(3,pp)
  xl(3,pp) = z_test

  do ll = 1, lint
    call FemShapeFunctions(sv(1,ll), xl, numDimensions, nel, xsj, shp)

    xsj = xsj*sv(5,ll)

    vol_test = vol_test + xsj
  enddo

  if ((vol_test<0.0d0).OR.(vol_test>volel)) then
    !write(6,*) "Point lies outside the tetrahedron"
    exit
  else
    kk = kk+1
  endif

  ! Give back xl nodes their original coordinates
  xl(1,pp) = temp_x
  xl(2,pp) = temp_y
  xl(3,pp) = temp_z
enddo

if (kk==4) then
  !write(6,*) "Point lies inside the tetrahedron"
  inside = 1
endif
!--------------------------------------------------------------------------------------------------------------------------------------------!
end subroutine FemIsPointInElement
