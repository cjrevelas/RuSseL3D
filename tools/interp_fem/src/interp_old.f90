subroutine interp(numnp, numel, num_plot_points, coef_x, coef_y, coef_z, coef_r, x_plot, y_plot, z_plot, r_plot, &
                                                                   & delta_plot, ix, xc, phi_matrix, phi_grafted, wa)
!--------------------------------------------------------------------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------------------------------------------------------------------!
integer, intent(in)                      :: numnp, numel, num_plot_points
integer, intent(in), dimension(4, numel) :: ix
integer                                  :: i, ii, kk, j, l, lint, n, m1, m2, inside

real(8), intent(in), dimension(3, numnp) :: xc
real(8), intent(in), dimension(numnp)    :: phi_matrix, phi_grafted, wa
real(8), intent(in)                      :: coef_x, coef_y, coef_z, coef_r
real(8), intent(inout)                   :: x_plot, y_plot, z_plot, r_plot, delta_plot
real(8), dimension(4,11)                 :: shp
real(8), dimension(5,11)                 :: sv
real(8), dimension(3,4)                  :: xl
real(8), dimension(4)                    :: ul_matrix, ul_grafted, ul_wa
real(8), dimension(4)                    :: gc, lc
real(8), dimension(4,4)                  :: transf, transf_inv
real(8)                                  :: volel, xsj
real(8)                                  :: u_interp_matrix, u_interp_grafted, u_interp_wa
!--------------------------------------------------------------------------------------------------------------------------------------------!
interface
    pure function matinv4(transf)
    real(8), intent(in), dimension(4,4) :: transf
    real(8),             dimension(4,4) :: matinv4
    end function
end interface
!--------------------------------------------------------------------------------------------------------------------------------------------!
open(unit = 80, file = "o.interp")

!set up for 11-point quadrature
l = 3
call gausspoints(l, lint, sv)

!loop over plot points
do kk = 1, num_plot_points

    x_plot = x_plot - (delta_plot / dble(num_plot_points)) * coef_x
    y_plot = y_plot - (delta_plot / dble(num_plot_points)) * coef_y
    z_plot = z_plot - (delta_plot / dble(num_plot_points)) * coef_z

    !r_plot = r_plot + (delta_plot / dble(num_plot_points)) * coef_r
    !x_plot = r_plot * sin(0.955316618) * cos(0.785398163)
    !y_plot = r_plot * sin(0.955316618) * sin(0.785398163)
    !z_plot = r_plot * cos(0.955316618)

    !loop over elements
    do n = 1, numel

        !loop over all nodes of current element
        do i = 1, 4
            !find global index of node
            ii = ix(i,n)

            !copy coordinates from global array to local array concerning current element
            do j = 1, 3
                xl(j,i) = xc(j,ii)
            enddo

            !copy value of solution at current node from global to local array
            ul_matrix(i)  = phi_matrix(ii)
            ul_grafted(i) = phi_grafted(ii)
            ul_wa(i)      = wa(ii)
        enddo

        !initialize accumulator for element volume
        volel = 0.d00

        !loop over all quadrature points in the considered element
        do l = 1, lint
            call tetshp(sv(1,l), xl, xsj, shp)

            xsj = xsj*sv(5,l)

            volel = volel + xsj
        enddo
        !-------------------------------------------START INTERPOLATION--------------------------------------------------------!
        !check if current plotPoint lies in the specific element
        call isInside(xl, x_plot, y_plot, z_plot, lint, sv, volel, inside, shp)

        if (inside==1) then
            gc(1) = 1.d00
            gc(2) = x_plot
            gc(3) = y_plot
            gc(4) = z_plot

            do m1 = 1, 4
                transf(1,m1) = 1.d00
                transf(2,m1) = xl(1,m1)
                transf(3,m1) = xl(2,m1)
                transf(4,m1) = xl(3,m1)
            enddo

            transf_inv = matinv4(transf)

            do m1 = 1, 4
                lc(m1) = 0.d00
                do m2 = 1, 4
                    lc(m1) = lc(m1) + transf_inv(m1,m2) * gc(m2)
                enddo
            enddo

            call tetshp(lc, xl, xsj, shp)

            u_interp_matrix  = 0.d00
            u_interp_grafted = 0.d00
            u_interp_wa      = 0.d00

            do j = 1, 4
                u_interp_matrix  = u_interp_matrix  + shp(4,j) * ul_matrix(j)
                u_interp_grafted = u_interp_grafted + shp(4,j) * ul_grafted(j)
                u_interp_wa      = u_interp_wa      + shp(4,j) * ul_wa(j)
            enddo

            write(80,'(I8,6(2X,E16.9))') n, x_plot, y_plot, z_plot, u_interp_matrix, u_interp_grafted, u_interp_wa
            goto 300
        endif
        !end loop over elements
        enddo
!end loop over plot points
300 enddo

close(80)
!--------------------------------------------------------------------------------------------------------------------------------------------!
end subroutine interp
