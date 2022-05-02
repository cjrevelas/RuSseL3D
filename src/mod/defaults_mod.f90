module defaults_mod
!----------------------------------------------------------------------------------------------------------------------------!
character(len=15), parameter :: dflt_inputFile = "./in.input"
character(len=15), parameter :: dflt_meshFile  = "./in.mesh"
character(len=15), parameter :: dflt_graftFile = "./in.gnodes"
character(len=15), parameter :: dflt_fieldFile = "./in.field.bin"

real(8), parameter :: dflt_ads_distance             = 1.2d1
real(8), parameter :: dflt_numGraftedChainsTol      = 5.0d-3
real(8), parameter :: dflt_r_gpoint                 = 4.0d-1
real(8), parameter :: dflt_bin_thickness            = 5.0d-1
real(8), parameter :: dflt_frac                     = 1.0d0
real(8), parameter :: dflt_fieldTol                 = 5.0d-1
real(8), parameter :: dflt_freeEnergyTol            = 1.0d-6
real(8), parameter :: dflt_freeEnergyTolForDelta    = 1.0d-3
real(8), parameter :: dflt_polymer_hamaker_constant = 0.0d0
real(8), parameter :: dflt_wall_distance            = 5.0d0
real(8), parameter :: dflt_pres                     = 0.0d0

integer, parameter :: dflt_iterations                  = 1
integer, parameter :: dflt_init_iter                   = 0
integer, parameter :: dflt_mx_exist                    = 0
integer, parameter :: dflt_gr_exist                    = 0
integer, parameter :: dflt_prof_dim                    = 3
integer, parameter :: dflt_mumps_matrix_type           = 0
integer, parameter :: dflt_export_phi_gen_freq         = 1
integer, parameter :: dflt_export_phi_indiv_freq       = 0
integer, parameter :: dflt_export_field_freq           = 1
integer, parameter :: dflt_export_field_bin_freq       = 1
integer, parameter :: dflt_export_propagators_freq     = 0
integer, parameter :: dflt_export_brush_thickness_freq = 1
integer, parameter :: dflt_export_chains_per_area_freq = 0
integer, parameter :: dflt_export_ads_free_freq        = 0

logical, parameter :: dflt_domain_is_periodic = .false.
logical, parameter :: dflt_square_gradient    = .false.
!----------------------------------------------------------------------------------------------------------------------------!
end module defaults_mod
