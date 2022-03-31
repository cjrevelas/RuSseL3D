!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module constants_mod
!----------------------------------------------------------------------------------------------------------------------------!
real(8), parameter :: boltz_const_Joule_molK = 8.3144598d0
real(8), parameter :: boltz_const_Joule_K    = 1.38064852d-23
real(8), parameter :: n_avog                 = 6.02214085700d23
real(8), parameter :: pi                     = 3.14159265359d0
real(8), parameter :: atm_to_pa              = 1.01325d5
real(8), parameter :: tol                    = 1.0d-12
real(8), parameter :: A_to_m                 = 1.0d-10
real(8), parameter :: A2_to_m2               = 1.0d-20
real(8), parameter :: A3_to_m3               = 1.0d-30
real(8), parameter :: m_to_A                 = 1.0d10
real(8), parameter :: m3_to_A3               = 1.0d30
real(8), parameter :: m3_to_cm3              = 1.0d6
real(8), parameter :: gr_cm3_to_kg_m3        = 1.0d3
real(8), parameter :: kg_m3_to_gr_m3         = 1.0d3
real(8), parameter :: N_to_mN                = 1.0d3
!----------------------------------------------------------------------------------------------------------------------------!
end module constants_mod
