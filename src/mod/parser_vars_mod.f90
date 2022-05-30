!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module parser_vars_mod
!----------------------------------------------------------------------------------------------------------------------------!
! Chain contour discretization
integer :: numEdwPointsMatrix, numConvolPointsMatrix
integer :: numEdwPointsGrafted, numConvolPointsGrafted
integer :: contourMatrix, contourGrafted
real(8) :: lengthMatrix, lengthMatrixMax, lengthGrafted
real(8) :: stepEdwAveMatrix, stepConvolAveMatrix
real(8) :: stepEdwAveGrafted, stepConvolAveGrafted
real(8) :: xs_crit_mx, xs_crit_gr

! Indices and flags
integer :: iow, ioe, iterations, fieldInitScheme, initialIterationId
integer :: matrixExist, graftedExist, deltaUpdateFreq, getICfromDelta
integer :: exportPhiGeneral, exportPhiIndividual, exportField, exportPropagators
integer :: exportFieldBinary, exportBrushThickness, exportChainsPerArea, exportAdsorbedFree
logical :: squareGradient

! Domain
integer               :: profileDimensions
integer               :: periodicity
integer, dimension(6) :: periodicFaceId
logical               :: domainIsPeriodic
logical, dimension(3) :: periodicAxisId
real(8)               :: graftPointDistance, binThickness

! Dirichlet boundaries
integer                              :: numDirichletFaces, numNanoparticleFaces
integer, allocatable, dimension(:)   :: dirichletFaceId, nanoparticleFaceId
real(8), allocatable, dimension(:)   :: dirichletFaceValue, nanoparticleFaceValue
real(8), allocatable, dimension(:)   :: radius_np_eff, A_np, sigma_np
real(8), allocatable, dimension(:)   :: A_plate, sigma_plate
real(8), allocatable, dimension(:,:) :: center_np

! SCF model and potential data
real(8) :: temperature, beta, pressure, massOfMonomer, massDensity, kapa, kappa_T
real(8) :: massBulkDensity, molarBulkDensity, segmentBulkDensity
real(8) :: k_gr, k_gr_tilde
real(8) :: rg2OfMatrixMonomer, rg2OfGraftedMonomer, sphere_radius
real(8) :: A_pol, sigma_pol, wall_distance, ads_distance

! Convergence scheme
integer :: mumpsMatrixType
real(8) :: fieldTol, freeEnergyTol, freeEnergyTolForDelta
real(8) :: frac, numGraftedChainsTol
!----------------------------------------------------------------------------------------------------------------------------!
end module parser_vars_mod
