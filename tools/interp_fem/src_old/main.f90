!--------------------------------------------------------------------------------------------------------------------------------------------!
!                  THIS IS A POST-PROCESSING CODE TO INTERPOLATE THE PROFILES OBTAINED BY THE FEM3D CODE                                     !
!                                 IT APPLIES A SPATIAL THREE-DIMENSIONAL INTERPOLATION USING                                                 !
!                          THE SAME FINITE-ELEMENT SHAPE FUNCTIONS AND MESH WHICH ARE USED IN FEM3D.                                         !
!                                                    DATE: 01/03/2020                                                                        !
!                      AUTHORS: CONSTANTINOS J. REVELAS, ARISTOTELIS P. SGOUROS, APOSTOLIS T. LAKKAS                                         !
!--------------------------------------------------------------------------------------------------------------------------------------------!
program fem3d_spatial_interpolation
!--------------------------------------------------------------------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------------------------------------------------------------------!
integer                              :: i, j, ii, numnp, numel, num_plot_points, dummy
integer, allocatable, dimension(:,:) :: ix

logical :: plot_across_x, plot_across_y, plot_across_z, plot_across_r

real(8)                              :: delta_plot, coef_x, coef_y, coef_z, coef_r
real(8)                              :: x_plot, x_plot_min, x_plot_max
real(8)                              :: y_plot, y_plot_min, y_plot_max
real(8)                              :: z_plot, z_plot_min, z_plot_max
real(8)                              :: r_plot, r_plot_min, r_plot_max
real(8)                              :: time_1, time_2, delta_time
real(8), allocatable, dimension(:)   :: phi_matrix, phi_grafted, wa
real(8), allocatable, dimension(:,:) :: xc
!--------------------------------------------------------------------------------------------------------------------------------------------!
call cpu_time(time_1)

open(unit = 22, file = "in.input.txt")
read(22,*) numnp
read(22,*) numel
read(22,*) num_plot_points

read(22,*)

read(22,*) plot_across_x
read(22,*) plot_across_y
read(22,*) plot_across_z
read(22,*) plot_across_r

read(22,*)

read(22,*) x_plot_min, x_plot_max
read(22,*) y_plot_min, y_plot_max
read(22,*) z_plot_min, z_plot_max
read(22,*) r_plot_min, r_plot_max
close(22)

if (plot_across_x) then
    coef_x = 1.0
    coef_y = 0.0
    coef_z = 0.0
    coef_r = 0.0

    delta_plot = x_plot_max - x_plot_min
    x_plot     = x_plot_min
    y_plot     = y_plot_min
    z_plot     = z_plot_min
    r_plot     = r_plot_min
elseif (plot_across_y) then
    coef_x = 0.0
    coef_y = 1.0
    coef_z = 0.0
    coef_r = 0.0

    delta_plot = y_plot_max - y_plot_min
    x_plot     = x_plot_min
    y_plot     = y_plot_min
    z_plot     = z_plot_min
    r_plot     = r_plot_min
elseif (plot_across_z) then
    coef_x = 0.0
    coef_y = 0.0
    coef_z = 1.0
    coef_r = 0.0

    delta_plot = z_plot_max - z_plot_min
    x_plot     = x_plot_min
    y_plot     = y_plot_min
    z_plot     = z_plot_max
    r_plot     = r_plot_min
elseif (plot_across_r) then
    coef_x = 0.0
    coef_y = 0.0
    coef_z = 0.0
    coef_r = 1.0

    delta_plot = r_plot_max - r_plot_min
    x_plot     = x_plot_min
    y_plot     = y_plot_min
    z_plot     = z_plot_min
    r_plot     = r_plot_min
endif

allocate(xc(3,numnp), ix(4,numel), phi_matrix(numnp), phi_grafted(numnp), wa(numnp))

open (unit = 77, file = "in.elemcon.txt")
read (77,*) ((ix(i,ii), i=1, 4), ii=1, numel)
close(77)

ix = ix + 1

open(unit = 88, file = "in.profiles.txt")
do i = 1, numnp
    read(88,*) dummy, (xc(j,i), j=1,3), phi_matrix(i), phi_grafted(i), wa(i)
enddo
close(88)

call interp(numnp, numel, num_plot_points, coef_x, coef_y, coef_z, coef_r, x_plot, y_plot, z_plot, &
                                                    & r_plot, delta_plot, ix, xc, phi_matrix, phi_grafted, wa)

call cpu_time(time_2)

delta_time = time_2 - time_1

write(6, '("cpu time: ", F12.2, "min")') delta_time/6.d01

deallocate(xc, ix, phi_matrix, phi_grafted, wa)
!--------------------------------------------------------------------------------------------------------------------------------------------!
end program fem3d_spatial_interpolation
