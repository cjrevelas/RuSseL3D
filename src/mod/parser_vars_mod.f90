!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module parser_vars_mod
!----------------------------------------------------------------------------------------------------------------------------!
! Chain contour discretization
integer :: ns_mx_ed, ns_mx_conv, ns_gr_ed, ns_gr_conv
integer :: contour_discr_mx, contour_discr_gr
real(8) :: chainlen_mx, chainlen_mx_max, chainlen_gr
real(8) :: ds_ave_mx_ed, ds_ave_gr_ed, ds_ave_mx_conv, ds_ave_gr_conv
real(8) :: xs_crit_mx, xs_crit_gr

! Indices and flags
integer :: iow, ioe, iterations, field_init_scheme, init_iter
integer :: mx_exist, gr_exist, calc_delta_every, grafted_ic_from_delta
integer :: export_phi_gen_freq, export_phi_indiv_freq, export_field_freq, export_propagators_freq
integer :: export_field_bin_freq, export_brush_thickness_freq, export_chains_per_area_freq, export_ads_free_freq
logical :: square_gradient

! Domain
integer               :: prof_dim
integer               :: periodicity
integer, dimension(6) :: periodic_face_id
logical               :: domain_is_periodic
logical, dimension(3) :: periodic_axis_id
real(8)               :: r_gpoint, bin_thickness

! Dirichlet boundaries
integer                              :: num_of_dirichlet_faces, num_of_nanoparticle_faces
integer, allocatable, dimension(:)   :: ids_dirichlet_faces, ids_nanopart_faces
real(8), allocatable, dimension(:)   :: radius_np_eff, A_np, sigma_np
real(8), allocatable, dimension(:)   :: A_plate, sigma_plate
real(8), allocatable, dimension(:,:) :: center_np

! SCF model and potential data
real(8) :: Temp, beta, Pres, mon_mass, massden, kapa, kappa_T
real(8) :: rho_mass_bulk, rho_mol_bulk, rho_seg_bulk
real(8) :: k_gr, k_gr_tilde
real(8) :: Rg2_per_mon_mx, Rg2_per_mon_gr, sphere_radius
real(8) :: A_pol, sigma_pol, wall_distance, ads_distance

! Convergence scheme
integer :: mumps_matrix_type
real(8) :: fieldTol, freeEnergyTol, freeEnergyTolForDelta
real(8) :: frac, numGraftedChainsTol
!----------------------------------------------------------------------------------------------------------------------------!
end module parser_vars_mod
