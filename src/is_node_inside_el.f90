!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine is_node_inside_el(xl, x_test, y_test, z_test, lint, sv, volel, inside, shp, ndm, nel)
!--------------------------------------------------------------------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------------------------------------------------------------------!
integer, intent(in)  :: lint, ndm, nel
integer, intent(out) :: inside
integer              :: p, k, l

real(8), intent(in)                     :: x_test, y_test, z_test
real(8), intent(in), dimension(5, 11)   :: sv
real(8), intent(inout), dimension(3, 4) :: xl
real(8), intent(out), dimension(4, 11)  :: shp
real(8), intent(out)                    :: volel
real(8)                                 :: temp_x, temp_y, temp_z, vol_test, xsj
!--------------------------------------------------------------------------------------------------------------------------------------------!
! Replace each coordinate vertex with the test point and calculate the new volume each time
inside = 0
k      = 0

!open(unit = 60, file = "d.elem_vol")
do p = 1, 4
    vol_test = 0.0d0

    temp_x  = xl(1,p)
    xl(1,p) = x_test

    temp_y  = xl(2,p)
    xl(2,p) = y_test

    temp_z  = xl(3,p)
    xl(3,p) = z_test

    do l = 1, lint
        call fem_tetshpfun(sv(1,l), xl, ndm, nel, xsj, shp)

        xsj = xsj*sv(5,l)

        vol_test = vol_test + xsj
    enddo

    if ((vol_test<0.).OR.(vol_test>volel)) then
         !write(6,*) "Point lies outside the tetrahedron"
         exit
    else
         k = k+1
    endif

    ! Give back xl nodes their original coordinates
    xl(1,p) = temp_x
    xl(2,p) = temp_y
    xl(3,p) = temp_z
enddo

if (k==4) then
    !write(6,*) "Point lies inside the tetrahedron"
    inside = 1
endif
!--------------------------------------------------------------------------------------------------------------------------------------------!
end subroutine is_node_inside_el
