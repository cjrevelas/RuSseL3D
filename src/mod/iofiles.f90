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
character(len=30) :: A_matrix_full   = "A_matrix_full.out.txt"
character(len=30) :: energy_terms    = "energy_terms.out.txt"
character(len=20) :: time            = "time.out.txt"
character(len=20) :: usolid          = "usolid.out.txt"
character(len=20) :: logfile         = "log.out.txt"
character(len=20) :: errorfile       = "error.out.txt"
character(len=30) :: matrix_assembly = "matrix_assembly_kcw.out.txt"
character(len=20) :: dir_faces       = "dir_faces.out.txt"
character(len=20) :: com_12          = "com_12.out.txt"
character(len=20) :: inter           = "inter.out.txt"
character(len=20) :: mesh_out        = "mesh.out.txt"
character(len=20) :: profiles        = "profiles.out.txt"
character(len=30) :: quad_interp     = "quad_interp_coeffs.out.txt"
character(len=30) :: simpson         = "simpson_coeffs.out.txt"
character(len=30) :: delta_out       = "delta.out.txt"
!-----------------------------------------------------------------------------------------------------------!
end module iofiles
