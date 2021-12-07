module parser_vars
!-------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------!
!chain contour discretizatiob
integer :: ns_matrix_ed, ns_matrix_conv, ns_gr_ed, ns_gr_conv
integer :: contour_integration_scheme
real(8) :: chainlen_matrix, chainlen_matrix_max, chainlen_gr
real(8) :: ds_ave_matrix_ed, ds_ave_gr_ed, ds_ave_matrix_conv, ds_ave_gr_conv

!indices and flags
integer :: iow, ioe, iterations, field_init_scheme, init_iter
integer :: matrix_exist, grafted_exist, calc_delta_every, grafted_ic_from_delta
integer :: print_ev, output_every
logical :: square_gradient

!domain
integer               :: prof_dim
real(8)               :: interf_area
real(8), dimension(3) :: box_lo, box_hi, box_len

!dirichlet boundaries
integer                              :: n_dirichlet_faces, n_nanopart_faces
integer, allocatable, dimension(:)   :: ids_dirichlet_faces, ids_nanopart_faces
real(8), allocatable, dimension(:)   :: radius_np, A_np, sigma_np
real(8), allocatable, dimension(:)   :: A_plate, sigma_plate
real(8), allocatable, dimension(:,:) :: center_np

!scf model and potential data
real(8) :: Temp, Pres, mon_mass, massden, kapa, kappa_T
real(8) :: rho_mass_bulk, rho_mol_bulk, rho_seg_bulk
real(8) :: k_gr, k_gr_tilde
real(8) :: Rg2_per_mon_matrix, Rg2_per_mon_gr, sphere_radius
real(8) :: A_pol, sigma_pol, wall_distance

!convergence scheme
integer :: scheme_type
integer :: mumps_matrix_type
real(8) :: max_error_tol, frac

!-------------------------------------------------------------------------------!
end module parser_vars
