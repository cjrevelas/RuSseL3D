module parser_vars
!-------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------!
!****************************TIME DISCRETIZATION********************************!
integer :: ns_free, ns_gr, ns_free_max, time_integration_scheme
real(8) :: chainlen_free, chainlen_gr
real(8) :: ds_ave_free, ds_ave_gr

!*******************************INDICES AND FLAGS*******************************!
integer :: iow, ioe, iterations, field_init_scheme, init_iter, output_every, print_ev
integer :: use_grafted

!************************************DOMAIN*************************************!
integer               :: prof_dim
real(8)               :: volume, interf_area
real(8), dimension(3) :: box_lo, box_hi, box_len

!*******************************BOUNDARY CONDITIONS*****************************!
integer                              :: n_dirichlet_faces, n_nanopart_faces
integer, allocatable, dimension(:)   :: ids_dirichlet_faces, ids_nanopart_faces
real(8), allocatable, dimension(:)   :: radius_np, A_np, sigma_np
real(8), allocatable, dimension(:)   :: A_plate, sigma_plate
real(8), allocatable, dimension(:,:) :: center_np

!*************************SCF MODEL AND POTENTIAL DATA**************************!
real(8) :: Temp, mon_mass, massden, kapa, kappa_T, rho_0
real(8) :: Rg2_per_mon_free, Rg2_per_mon_gr, sphere_radius
real(8) :: A_pol, sigma_pol

!*************************SCF CONVERGENCE SCHEME********************************!
integer :: scheme_type
real(8) :: max_error_tol, frac, mix_coef_frac, mix_coef_kapa

!*********************************MUMPS OPTIONS*********************************!
integer :: mumps_matrix_type
!-------------------------------------------------------------------------------!
end module parser_vars
