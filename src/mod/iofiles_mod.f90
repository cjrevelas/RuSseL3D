!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module iofiles_mod
!----------------------------------------------------------------------------------------------------------------------------!
! Input files
character(len=50) :: IO_inputFile
character(len=50) :: IO_meshFile
character(len=50) :: IO_graftFile
character(len=50) :: IO_fieldFile

! Error files
character(len=20) :: IO_errorFile = "e.error"

! Debug files
character(len=30) :: IO_stiffnessMatrix  = "d.stiffness_matrix"
character(len=30) :: IO_contourCoeffs    = "d.contour_coeffs"
character(len=30) :: IO_matrixAssembly   = "d.matrix_assembly"
character(len=20) :: IO_dirichletFaces   = "d.dir_faces"
character(len=20) :: IO_nodePairs        = "d.com_12"
character(len=20) :: IO_nodeConnectivity = "d.inter"
character(len=20) :: IO_nodeCoordinates  = "d.mesh"
character(len=20) :: IO_meshProfile      = "d.prof_mp"
character(len=20) :: IO_vertexElements   = "d.vertex_elements"
character(len=20) :: IO_edgeElements     = "d.edge_elements"
character(len=20) :: IO_faceElements     = "d.face_elements"
character(len=20) :: IO_domainElements   = "d.domain_elements"
character(len=20) :: IO_xFace1Elements   = "d.xface1_elements"
character(len=20) :: IO_xFace2Elements   = "d.xface2_elements"
character(len=20) :: IO_yFace1Elements   = "d.yface1_elements"
character(len=20) :: IO_yFace2Elements   = "d.yface2_elements"
character(len=20) :: IO_zFace1Elements   = "d.zface1_elements"
character(len=20) :: IO_zFace2Elements   = "d.zface2_elements"
character(len=20) :: IO_nodePairingXX    = "d.node_pairing_xx"
character(len=20) :: IO_nodePairingYY    = "d.node_pairing_yy"
character(len=20) :: IO_nodePairingZZ    = "d.node_pairing_zz"

! Output files
character(len=20) :: IO_solidPotential = "o.usolid"
character(len=20) :: IO_nodalProfile   = "o.phi"
character(len=20) :: IO_VtuProfiles    = "o.profiles.vtu"
character(len=20) :: IO_VtuGrafted     = "o.grafted.vtu"
character(len=20) :: IO_radialProfile  = "o.phi_radial"
character(len=20) :: IO_indivProfile   = "o.phi_indiv"
character(len=20) :: IO_field          = "o.field"
character(len=30) :: IO_delta          = "o.delta"
character(len=30) :: IO_brush          = "o.brush"
character(len=30) :: IO_energies       = "o.energy_terms"
character(len=20) :: IO_logFile        = "o.log"
character(len=20) :: IO_time           = "o.time"
!----------------------------------------------------------------------------------------------------------------------------!
end module iofiles_mod
