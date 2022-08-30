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
real(8) :: critContourMatrix, critContourGrafted

! Indices and flags
logical :: squareGradient
integer :: iow, ioe, iterations, fieldInitScheme, initialIterationId
integer :: matrixExist, graftedExist, deltaUpdateFreq, getICfromDelta
integer :: exportPhiGeneral, exportPhiIndividual, exportField, exportPropagators
integer :: exportFieldBin, exportBrushThickness, exportChainsPerArea, exportAdsorbedFree
integer :: exportAllGraftedChains
integer :: numGraftedChainsToExport
integer, allocatable, dimension(:) :: gpIndexToExport

! Domain
integer               :: profileDimensions
integer               :: periodicity
integer, dimension(6) :: periodicFaceId
logical               :: domainIsPeriodic
logical, dimension(3) :: periodicAxisId
real(8)               :: graftPointDistance, binThickness

! Dirichlet boundaries
integer                              :: numDirichletFaces, numNanopFaces
integer, allocatable, dimension(:)   :: dirichletFaceId, nanopFaceId
real(8), allocatable, dimension(:)   :: dirichletFaceValue, nanopFaceValue
real(8), allocatable, dimension(:)   :: nanopRadiusEff, nanopAlpha, nanopSigma
real(8), allocatable, dimension(:)   :: plateAlpha, plateSigma
real(8), allocatable, dimension(:,:) :: nanopCenter

! SCF model and potential data
real(8) :: temperature, beta, pressure, massOfMonomer, massDensity, kapa
real(8) :: massBulkDensity, molarBulkDensity, segmentBulkDensity
real(8) :: sgtParam, sgtParamTilde
real(8) :: rg2OfMatrixMonomer, rg2OfGraftedMonomer
real(8) :: polymerAlpha, polymerSigma
real(8) :: wallDistance, adsorptionDistance

! Convergence scheme
integer :: mumpsMatrixType
real(8) :: fieldTol, freeEnergyTol, freeEnergyTolForDelta
real(8) :: frac, numGraftedChainsTol
!----------------------------------------------------------------------------------------------------------------------------!
end module parser_vars_mod
