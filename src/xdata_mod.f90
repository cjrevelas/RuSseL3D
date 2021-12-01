module xdata
!-------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------!
!****************************TIME DISCRETIZATION********************************!
integer                            :: ns, time_integration_scheme
real(8)                            :: chainlen
real(8)                            :: ds_ave

!*******************************INDICES AND FLAGS*******************************!
integer :: iow, ioe, iterations, readfield, init_iter, output_every, print_ev
integer :: geom_type, use_grafted

!************************************DOMAIN*************************************!
integer               :: prof_dim
real(8)               :: volume
real(8), dimension(3) :: box_lo, box_hi, box_len

!*******************************BOUNDARY CONDITIONS*****************************!
integer                            :: n_dirichlet_faces
integer, allocatable, dimension(:) :: ids_dirichlet_faces

!*************************SCF MODEL AND POTENTIAL DATA**************************!
real(8) :: Temp, mon_mass, massden, kapa, kappa_T, rho_0, Rg_2, CN
real(8) :: Aps, Asio2, sigma1, sigma2, sphere_radius, distance

!*************************SCF CONVERGENCE SCHEME********************************!
real(8) :: max_error_tol, max_error, std_error, frac, mix_coef_frac, mix_coef_kapa
integer :: scheme_type

!*********************************MUMPS OPTIONS*********************************!
integer :: mumps_matrix_type
!-------------------------------------------------------------------------------!
end module xdata
