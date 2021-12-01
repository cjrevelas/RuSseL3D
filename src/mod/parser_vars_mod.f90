module parser_vars
!-------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------!
!****************************TIME DISCRETIZATION********************************!
integer :: ns_free, ns_gr, ns_free_max, time_integration_scheme
real(8) :: chainlen_free, chainlen_gr
real(8) :: ds_ave_free, ds_ave_gr

!*******************************INDICES AND FLAGS*******************************!
integer :: iow, ioe, iterations, readfield, init_iter, output_every, print_ev
integer :: geom_type, use_grafted

!************************************DOMAIN*************************************!
integer               :: prof_dim
real(8)               :: volume, interf_area                                                         !CJR
real(8), dimension(3) :: box_lo, box_hi, box_len

!*******************************BOUNDARY CONDITIONS*****************************!
integer                            :: n_dirichlet_faces
integer, allocatable, dimension(:) :: ids_dirichlet_faces

!*************************SCF MODEL AND POTENTIAL DATA**************************!
real(8) :: Temp, mon_mass, massden, kapa, kappa_T, rho_0
real(8) :: Rg2_per_mon_free, Rg2_per_mon_gr, sphere_radius
real(8) :: Aps, Asio2, sigma1, sigma2

!*************************SCF CONVERGENCE SCHEME********************************!
integer :: scheme_type
real(8) :: max_error_tol, frac, mix_coef_frac, mix_coef_kapa

!*********************************MUMPS OPTIONS*********************************!
integer :: mumps_matrix_type
!-------------------------------------------------------------------------------!
end module parser_vars
