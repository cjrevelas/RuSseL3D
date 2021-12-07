module constants
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
real(8), parameter :: boltz_const_Joule_molK = 8.3144598         !Boltzmann constant in J/mol/K.
real(8), parameter :: boltz_const_Joule_K    = 1.38064852E-23    !Boltzmann constant in J/K.
real(8), parameter :: n_avog                 = 6.02214085700E+23
real(8), parameter :: pi                     = 3.14159265359
real(8), parameter :: tol                    = 1.E-12
real(8), parameter :: gr_cm3_to_kg_m3        = 1000              !kg/m^3
real(8), parameter :: atm_to_pa              = 1.01325d05        !pa
!--------------------------------------------------------------------!
end module constants
