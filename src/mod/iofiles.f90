!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module iofiles
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
!input files
character(len=50) :: input_filename
character(len=50) :: mesh_filename
character(len=50) :: gp_filename
character(len=50) :: field_in_filename

!error files
character(len=20) :: errorfile = "e.error"

!debug files
character(len=30) :: A_matrix_full   = "d.stiffness_matrix"
character(len=30) :: contour_coeffs  = "d.contour_coeffs"
character(len=30) :: matrix_assembly = "d.matrix_assembly"
character(len=20) :: dir_faces       = "d.dir_faces"
character(len=20) :: com_12          = "d.com_12"
character(len=20) :: inter           = "d.inter"
character(len=20) :: mesh_out        = "d.mesh"
character(len=20) :: mesh_prof       = "d.prof_mp"

!output files
character(len=20) :: usolid            = "o.usolid"
character(len=20) :: phi_nodal         = "o.phi"
character(len=20) :: phi_radial        = "o.phi_radial"
character(len=20) :: phi_profile_indiv = "o.phi_indiv"
character(len=20) :: field_profile     = "o.field"
character(len=30) :: delta_out         = "o.delta"
character(len=30) :: brush_out         = "o.brush"
character(len=30) :: energy_terms      = "o.energy_terms"
character(len=20) :: logfile           = "o.log"
character(len=20) :: time              = "o.time"
!-----------------------------------------------------------------------------------------------------------!
end module iofiles
