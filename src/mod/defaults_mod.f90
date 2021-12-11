module defaults_mod
    character(len=15), parameter :: default_input_filename    = "./in.input"
    character(len=15), parameter :: default_mesh_filename     = "./in.mesh"
    character(len=15), parameter :: default_gp_filename       = "./in.gnodes"
    character(len=15), parameter :: default_field_in_filename = "./in.field.bin"

    real(8), parameter :: default_ads_distance                    = 1.2d1
    real(8), parameter :: default_num_gr_chains_tol               = 5.0d-3
    real(8), parameter :: default_r_gpoint                        = 4.0d-1
    real(8), parameter :: default_bin_thickness                   = 5.0d-1
    real(8), parameter :: default_frac                            = 1.0d0
    real(8), parameter :: default_field_error_tol                 = 5.0d-1
    real(8), parameter :: default_free_energy_error_tol           = 1.0d-6
    real(8), parameter :: default_free_energy_error_tol_for_delta = 1.0d-3
    real(8), parameter :: default_polymer_hamaker_constant        = 0.0d0
    real(8), parameter :: default_wall_distance                   = 5.0d0
    real(8), parameter :: default_pres                            = 0.0d0

    integer, parameter :: default_iterations                  = 1
    integer, parameter :: default_init_iter                   = 0
    integer, parameter :: default_mx_exist                    = 0
    integer, parameter :: default_gr_exist                    = 0
    integer, parameter :: default_prof_dim                    = 3
    integer, parameter :: default_mumps_matrix_type           = 0
    integer, parameter :: default_export_phi_gen_freq         = 1
    integer, parameter :: default_export_phi_indiv_freq       = 0
    integer, parameter :: default_export_field_freq           = 1
    integer, parameter :: default_export_field_bin_freq       = 1
    integer, parameter :: default_export_propagators_freq     = 0
    integer, parameter :: default_export_brush_thickness_freq = 1
    integer, parameter :: default_export_chains_per_area_freq = 0
    integer, parameter :: default_export_ads_free_freq        = 0

    logical, parameter :: default_domain_is_periodic = .false.
    logical, parameter :: default_square_gradient    = .false.
end module defaults_mod
