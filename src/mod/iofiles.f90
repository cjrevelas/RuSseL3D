module iofiles
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
!input files
character(len=20) :: input_filename    = "input.in.txt"
character(len=20) :: mesh_filename     = "mesh.in.mphtxt"
character(len=20) :: gp_filename       = "gnodes.in.lammpstrj"
character(len=20) :: field_in_filename = "field.in.bin"

!output files
character(len=30) :: A_matrix_full      = "d.A_matrix_full.out.txt"
character(len=20) :: time               = "d.time.out.txt"
character(len=30) :: quad_interp        = "d.quad_interp_coeffs.out.txt"
character(len=30) :: simpson            = "d.simpson_coeffs.out.txt"
character(len=20) :: errorfile          = "d.error.out.txt"
character(len=30) :: matrix_assembly    = "d.matrix_assembly_kcw.out.txt"
character(len=20) :: dir_faces          = "d.dir_faces.out.txt"
character(len=20) :: com_12             = "d.com_12.out.txt"
character(len=20) :: inter              = "d.inter.out.txt"
character(len=20) :: mesh_out           = "d.mesh.out.txt"
character(len=20) :: mesh_prof          = "d.prof_mp.out.txt"
character(len=20) :: usolid             = "d.usolid.out.txt"
character(len=20) :: phi_nodal          = "o.phi.out.txt"
character(len=20) :: phi_radial         = "o.phi_radial.out.txt"
character(len=20) :: phi_profile_indiv  = "o.phi_indiv.out.txt"
character(len=20) :: field_profile      = "o.field.out.txt"
character(len=30) :: delta_out          = "o.delta.out.txt"
character(len=30) :: brush_out          = "o.brush.out.txt"
character(len=30) :: energy_terms       = "o.energy_terms.out.txt"
character(len=20) :: logfile            = "o.log.out.txt"
!-----------------------------------------------------------------------------------------------------------!
end module iofiles
