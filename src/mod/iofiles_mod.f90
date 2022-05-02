!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module iofiles_mod
!----------------------------------------------------------------------------------------------------------------------------!
! Input files
character(len=50) :: inputFile
character(len=50) :: meshFile
character(len=50) :: graftFile
character(len=50) :: fieldFile

! Error files
character(len=20) :: errorfile = "e.error"

! Debug files
character(len=30) :: A_matrix_full   = "d.stiffness_matrix"
character(len=30) :: contour_coeffs  = "d.contour_coeffs"
character(len=30) :: matrix_assembly = "d.matrix_assembly"
character(len=20) :: dir_faces       = "d.dir_faces"
character(len=20) :: com_12          = "d.com_12"
character(len=20) :: inter           = "d.inter"
character(len=20) :: mesh_out        = "d.mesh"
character(len=20) :: mesh_prof       = "d.prof_mp"
character(len=20) :: vertex_elements = "d.vertex_elements"
character(len=20) :: edge_elements   = "d.edge_elements"
character(len=20) :: face_elements   = "d.face_elements"
character(len=20) :: domain_elements = "d.domain_elements"
character(len=20) :: xface1_elements = "d.xface1_elements"
character(len=20) :: xface2_elements = "d.xface2_elements"
character(len=20) :: yface1_elements = "d.yface1_elements"
character(len=20) :: yface2_elements = "d.yface2_elements"
character(len=20) :: zface1_elements = "d.zface1_elements"
character(len=20) :: zface2_elements = "d.zface2_elements"
character(len=20) :: node_pairing_xx = "d.node_pairing_xx"
character(len=20) :: node_pairing_yy = "d.node_pairing_yy"
character(len=20) :: node_pairing_zz = "d.node_pairing_zz"

! Output files
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
!----------------------------------------------------------------------------------------------------------------------------!
end module iofiles_mod
