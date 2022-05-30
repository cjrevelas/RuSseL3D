module defaults_mod
!----------------------------------------------------------------------------------------------------------------------------!
use flags_mod, only: mumps_asymm
!----------------------------------------------------------------------------------------------------------------------------!
character(len=15), parameter :: dflt_inputFile = "./in.input"
character(len=15), parameter :: dflt_meshFile  = "./in.mesh"
character(len=15), parameter :: dflt_graftFile = "./in.gnodes"
character(len=15), parameter :: dflt_fieldFile = "./in.field.bin"

real(8), parameter :: dflt_ads_distance             = 1.2d1
real(8), parameter :: dflt_numGraftedChainsTol      = 5.0d-3
real(8), parameter :: dflt_graftPointDistance       = 4.0d-1
real(8), parameter :: dflt_binThickness             = 5.0d-1
real(8), parameter :: dflt_fraction                 = 1.0d0
real(8), parameter :: dflt_fieldTol                 = 5.0d-1
real(8), parameter :: dflt_freeEnergyTol            = 1.0d-6
real(8), parameter :: dflt_freeEnergyTolForDelta    = 1.0d-3
real(8), parameter :: dflt_polymer_hamaker_constant = 0.0d0
real(8), parameter :: dflt_wall_distance            = 5.0d0
real(8), parameter :: dflt_pressure                 = 0.0d0
real(8), parameter :: dflt_dirichletValue           = 0.0d0
real(8), parameter :: dflt_nanoparticleValue        = 0.0d0

integer, parameter :: dflt_iterations           = 1
integer, parameter :: dflt_initialIterationId   = 0
integer, parameter :: dflt_matrixExist          = 0
integer, parameter :: dflt_graftedExist         = 0
integer, parameter :: dflt_profileDimensions    = 3
integer, parameter :: dflt_mumpsMatrixType      = mumps_asymm
integer, parameter :: dflt_exportPhiGeneral     = 1
integer, parameter :: dflt_exportPhiIndividual  = 0
integer, parameter :: dflt_exportField          = 1
integer, parameter :: dflt_exportFieldBinary    = 1
integer, parameter :: dflt_exportPropagators    = 0
integer, parameter :: dflt_exportBrushThickness = 1
integer, parameter :: dflt_exportChainsPerArea  = 0
integer, parameter :: dflt_exportAdsorbedFree   = 0

logical, parameter :: dflt_domainIsPeriodic = .false.
logical, parameter :: dflt_squareGradient   = .false.
!----------------------------------------------------------------------------------------------------------------------------!
end module defaults_mod
