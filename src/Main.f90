!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

program RuSseL
!----------------------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod
use arrays_mod
use eos_mod
use constants_mod
use error_handling_mod
use write_helper_mod
use geometry_mod
use iofiles_mod
use delta_mod
use hist_mod
use kcw_mod, only: rdiag1, F_m
use flags_mod, only: contourSymmetric, contourAsymmetric, contourHybrid
use fhash_module__ints_double
use ints_module
#ifdef USE_MPI
use mpistuff_mod
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
#ifdef USE_MPI
include "mpif.h"
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
integer :: iter, ToolsSystemTime, tInit, tFinal, node, graftNodeId, graftedChain

logical :: convergence = .False., computeDelta = .False.

real(8), allocatable, dimension(:) :: dphi2_dr2

type(fhash_type__ints_double) :: elemcon

real(8) :: partitionMatrixChains = 0.0d0, numMatrixChains = 0.0d0
real(8) :: numGraftedChains = 0.0d0, numGraftedChainsError = 1.0d2
real(8) :: wwFieldMaximum = 0.0d0, wwFieldStdError = 0.0d0, wwFieldError = 2.0d5
real(8) :: freeEnergyPrevious = 1.0d10, freeEnergy = 0.0d0, freeEnergyError = 1.0d10
!----------------------------------------------------------------------------------------------------------------------------------!
!**************************************************************************************************************!
!                                                    MPI SECTION                                               !
!**************************************************************************************************************!
#ifdef USE_MPI
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

if (my_id==0) then
  root = .true.
else
  root = .False.
endif

if (root) then
  write (6,*)
  write (6,'("MPI run with ",I4," procs")') n_proc
  write (6,*)
endif

flag_continue = .true.

! The slave processes will enter the mumps subroutine until they receive a stop signal from master proc
if (.not.root) then
  ! Receive the matrix type from root
  CALL MPI_BCAST(mumpsMatrixType, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
  do while (.true.)
    CALL MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (flag_continue) then
      CALL SolverMumps(mumpsMatrixType)
    else
      exit
    endif
  end do
  goto 1000
endif
#endif

iow = 10
open(unit=iow, file = IO_logFile)

! Initialize the error log
ioe = 11
open(unit=ioe, file = IO_errorFile, status='replace')
close(ioe)
!**************************************************************************************************************!
!                                             INITIALIZATION SECTION                                           !
!**************************************************************************************************************!
CALL ParserInput()
CALL ParserMesh(elemcon)
CALL InitArrays()

do node = 1, numNodes
  CALL ComputeNodeVolume(nodeVolume(node), node)
enddo

allocate(dphi2_dr2(numNodes))
dphi2_dr2=0.0d0

CALL InitDelta()
CALL ToolsHistogram(binThickness, nodeVolume)

#ifdef USE_MPI
CALL MPI_BCAST(mumpsMatrixType, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif

CALL InitChainContour(contourMatrix, lengthMatrixMax, critContourMatrix, numEdwPointsMatrix, stepEdwAveMatrix, dsEdwMatrix, xsEdwMatrix, coeffEdwMatrix)

! In absence of mx chains, skip discretization for mx convolution
if (matrixExist.eq.1) then
  if (contourMatrix.ne.contourUniform) then
    CALL InitChainContour(contourSymmetric, lengthMatrix, critContourMatrix, numConvolPointsMatrix, stepConvolAveMatrix, dsConvMatrix, xsConvMatrix, coeffConvMatrix)
  else
    CALL InitChainContour(contourMatrix, lengthMatrix, critContourMatrix, numConvolPointsMatrix, stepConvolAveMatrix, dsConvMatrix, xsConvMatrix, coeffConvMatrix)
  endif
endif

if (graftedExist.eq.1) then
  CALL InitChainContour(contourGrafted, lengthGrafted, critContourGrafted, numEdwPointsGrafted, stepEdwAveGrafted, dsEdwGrafted, xsEdwGrafted, coeffEdwGrafted)
  if (contourGrafted.ne.contourUniform) then
    CALL InitChainContour(contourSymmetric, lengthGrafted, critContourGrafted, numConvolPointsGrafted, stepConvolAveGrafted, dsConvGrafted, xsConvGrafted, coeffConvGrafted)
  else
    CALL InitChainContour(contourGrafted, lengthGrafted, critContourGrafted, numConvolPointsGrafted, stepConvolAveGrafted, dsConvGrafted, xsConvGrafted, coeffConvGrafted)
  endif
endif

CALL InitField(uuField, wwField)

wwFieldMixed = wwField
!**************************************************************************************************************!
!                                        LOOPS FOR FIELD CONVERGENCE                                           !
!**************************************************************************************************************!
write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-----------------------------------SIMULATION STARTS-----------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------SIMULATION STARTS-----------------------------------',85)

write(iow,'(A10,1X,8(A19,1X),A16)') "iter", "fraction", "energy", "energy_error", "gr_chains", "gr_chains_error", "field_error", "field_std_error", "field_max"
write(6  ,'(A4,1X,8(A14,1X),A12)')  "iter", "fraction", "energy", "energy_error", "gr_chains", "gr_chains_error", "field_error", "field_std_error", "field_max"

tInit = ToolsSystemTime()

do iter = initialIterationId, iterations-1
  write(iow,'(I10,1X,8(E19.9E3,1X))') iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, wwFieldError, wwFieldStdError, wwFieldMaximum
  write(6  ,'(I4 ,1X,8(E14.4E3,1X))') iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, wwFieldError, wwFieldStdError, wwFieldMaximum

  close(iow)
  open(unit=iow, file = IO_logFile, position = 'append')

  wwField = wwFieldMixed

  CALL FemMatrixAssemble(rg2OfMatrixMonomer, wwField)

  do node = 1, numNodes
    qqMatrix(1,node)      = 1.0d0
    qqMatrixFinal(1,node) = 1.0d0
  enddo

  ! We always need to solve for mx chains, so that we have the ability to perform convolution
  CALL SolverEdwards(dsEdwMatrix, numEdwPointsMatrix, qqMatrix, qqMatrixFinal, nodeBelongsToDirichletFace, elemcon)

  if (graftedExist.eq.1) then
    do node = 1, numNodes
      CALL interp_linear(1, numEdwPointsMatrix+1, xsEdwMatrix, qqMatrixFinal(:,node), numConvolPointsGrafted+1, xsConvGrafted, qqMatrixInterpGrafted(:,node))
    enddo

    ! Recompute Delta functions, when the error in the number of grafted chains exceeds the specified tolerance
    if (getICfromDelta.eq.1) then
      computeDelta = ((iter==0) .or. ((freeEnergyError <= freeEnergyTolForDelta) .and. (numGraftedChainsError > numGraftedChainsTol)))

      if (computeDelta) then
        CALL ComputeDeltaNumerical(elemcon, qqMatrixInterpGrafted, wwFieldMixed, deltaNumerical)

        CALL ExportDelta(qqMatrixInterpGrafted, deltaNumerical, graftPointValue)
      endif

      do graftedChain = 1, targetNumGraftedChains
        graftNodeId = graftPointId(graftedChain)

        graftPointValue(graftedChain) = deltaNumerical(graftedChain) * lengthGrafted * 1.0d0 / &
                                        (qqMatrixInterpGrafted(numConvolPointsGrafted+1,graftNodeId) * (molarBulkDensity * n_avog))
      enddo
    endif

    CALL FemMatrixAssemble(rg2OfGraftedMonomer, wwField)

    qqGrafted      = 0.0d0
    qqGraftedFinal = 0.0d0

    do graftedChain = 1, targetNumGraftedChains
      graftNodeId = graftPointId(graftedChain)

      qqGrafted(1,graftNodeId)      = graftPointValue(graftedChain)
      qqGraftedFinal(1,graftNodeId) = graftPointValue(graftedChain)
    enddo

    CALL SolverEdwards(dsEdwGrafted, numEdwPointsGrafted, qqGrafted, qqGraftedFinal, nodeBelongsToDirichletFace, elemcon)
  endif

  ! In absence of mx chains, skip mx convolution
  if (matrixExist.eq.1) then
    do node = 1, numNodes
      CALL interp_linear(1, numEdwPointsMatrix+1, xsEdwMatrix, qqMatrixFinal(:,node), numConvolPointsMatrix+1, xsConvMatrix, qqMatrixInterp(:,node))
    enddo

    CALL ContourConvolution(lengthMatrix, numConvolPointsMatrix, coeffConvMatrix, qqMatrixInterp, qqMatrixInterp, phiMatrix)
  endif

  if (graftedExist.eq.1) then
    do node = 1, numNodes
      CALL interp_linear(1, numEdwPointsGrafted+1, xsEdwGrafted, qqGraftedFinal(:,node), numConvolPointsGrafted+1, xsConvGrafted, qqGraftedInterp(:,node))
    enddo

    CALL ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeffConvGrafted, qqGraftedInterp, qqMatrixInterpGrafted, phiGrafted)
  endif

  phiTotal = 0.0d0
  do node = 1, numNodes
    if (matrixExist.eq.1)  phiTotal(node) = phiTotal(node) + phiMatrix(node)
    if (graftedExist.eq.1) phiTotal(node) = phiTotal(node) + phiGrafted(node)
  enddo

  if (matrixExist.eq.1)  CALL ComputePartitionMatrix(numNodes, numConvolPointsMatrix, qqMatrixInterp, partitionMatrixChains)
  if (matrixExist.eq.1)  CALL ComputeNumberOfChains(lengthMatrix, phiMatrix, numMatrixChains)
  if (graftedExist.eq.1) CALL ComputeNumberOfChains(lengthGrafted, phiGrafted, numGraftedChains)

  ! Compare this to russel1d, to be sure about taking into account the bulk contribution or not
  do node = 1, numNodes
    wwFieldNew(node) = (eos_df_drho(phiTotal(node)) - eos_df_drho(1.0d0)) / (boltz_const_Joule_K*temperature) - &
                     &  sgtParam * (segmentBulkDensity * dphi2_dr2(node)) / (boltz_const_Joule_K * temperature) + uuField(node)
  enddo

  wwFieldError    = 0.0d0
  wwFieldStdError = 0.0d0
  wwFieldMaximum  = 0.0d0

  do node = 1, numNodes
    wwFieldError    = MAX(wwFieldError,DABS(wwFieldNew(node) - wwField(node)))
    wwFieldStdError = wwFieldStdError + (wwFieldNew(node) - wwField(node))**2.0d0
    wwFieldMaximum  = MAX(wwFieldMaximum, wwFieldNew(node))
  enddo

  wwFieldStdError = SQRT(wwFieldStdError / FLOAT((numNodes - 1)))
  wwFieldMaximum  = wwFieldMaximum  * lengthMatrix
  wwFieldError    = wwFieldError    * lengthMatrix
  wwFieldStdError = wwFieldStdError * lengthMatrix

  freeEnergyError    = ABS(freeEnergy - freeEnergyPrevious)
  freeEnergyPrevious = freeEnergy

  if (graftedExist.eq.1) then
    numGraftedChainsError = ABS(numGraftedChains-DBLE(targetNumGraftedChains)) / DBLE(targetNumGraftedChains)
  else
    numGraftedChainsError = 0.0d0
  endif

  do node = 1, numNodes
    wwFieldMixed(node) = (1.0d0 - frac) * wwField(node) + frac * wwFieldNew(node)
  enddo

  convergence = (wwFieldError<=fieldTol).or.(freeEnergyError<=freeEnergyTol)

  CALL ExportFieldBinary(wwFieldMixed, numNodes, 0)

  if (export(exportFieldBin, iter, convergence)) CALL ExportFieldBinary(wwFieldMixed, numNodes, iter)

  if ((MOD(iter,1).eq.0).OR.convergence) CALL ExportEnergies(qqMatrixInterpGrafted, qqGraftedInterp, phiTotal, wwFieldNew, uuField, partitionMatrixChains, targetNumGraftedChains, graftPointId, freeEnergy)

  CALL ExportComputes(iter, convergence, elemcon)

  CALL ExportVtuProfiles(phiMatrix, phiGrafted, wwFieldMixed)

  if (convergence) exit
enddo
!**************************************************************************************************************!
!                                             EXPORT SIMULATION RESULTS                                        !
!**************************************************************************************************************!
write(iow,'(I10,1X,8(E19.9E3,1X))')  iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, wwFieldError, wwFieldStdError, wwFieldMaximum
write(6  ,'(I4 ,1X,8(E14.4E3,1X))')  iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, wwFieldError, wwFieldStdError, wwFieldMaximum

write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-----------------------------------SUMMARIZED RESULTS-----------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------SUMMARIZED RESULTS-----------------------------------',85)

if (wwFieldError.lt.fieldTol) then
  write(iow,'("Field convergence of max error",F16.9)') wwFieldError
  write(6  ,'("Field convergence of max error",F16.9)') wwFieldError
endif

if (freeEnergyError.lt.freeEnergyTol) then
  write(iow,'("Energy convergence of max error",F16.9)') freeEnergyError
  write(6  ,'("Energy convergence of max error",F16.9)') freeEnergyError
endif

write(iow,'(3X,A40,E16.9)')adjl("Free energy (mJ/m2):",40), freeEnergy
write(6  ,'(3X,A40,E16.9)')adjl("Free energy (mJ/m2):",40), freeEnergy
write(iow,'(3X,A40,E16.9)')adjl("Interface area (A2):",40), interfaceArea()
write(6  ,'(3X,A40,E16.9)')adjl("Interface area (A2):",40), interfaceArea()
if (matrixExist.eq.1) then
  write(iow,'(3X,A40,E16.9)')adjl("Partition function of matrix chains:",40), partitionMatrixChains
  write(6  ,'(3X,A40,E16.9)')adjl("Partition function of matrix chains:",40), partitionMatrixChains
  write(iow,'(3X,A40,E16.4)')adjl("Number of matrix chains:",40),             numMatrixChains
  write(6  ,'(3X,A40,E16.4)')adjl("Number of matrix chains:",40),             numMatrixChains
endif
if (graftedExist.eq.1) then
  write(iow,'(3X,A40,E16.9)')adjl("Grafting density (A^-2):",40),  numGraftedChains/interfaceArea()
  write(6  ,'(3X,A40,E16.9)')adjl("Grafting density (A^-2):",40),  numGraftedChains/interfaceArea()
  write(iow,'(3X,A40,E16.4)')adjl("Number of grafted chains:",40), numGraftedChains
  write(6  ,'(3X,A40,E16.4)')adjl("Number of grafted chains:",40), numGraftedChains
endif

tFinal = ToolsSystemTime()
write(6,'(3X,A40,I16)')adjl('Run duration:',40), tFinal - tInit

! Deallocate all remaining dynamic memory
deallocate(nodeCoord)
deallocate(dphi2_dr2, d2phi_dr2)
if (numDirichletFaces > 0) deallocate(dirichletFaceId, dirichletFaceValue, plateAlpha, plateSigma)
if (numNanopFaces > 0) deallocate(nanopFaceId, nanopFaceValue, nanopAlpha, nanopSigma, nanopRadiusEff, nanopCenter)
if (periodicAxisId(1).and.periodicAxisId(2)) deallocate(edgeNodeOneXY, edgeNodeTwoXY, edgeNodeThreeXY, edgeNodeFourXY)
if (periodicAxisId(1).and.periodicAxisId(3)) deallocate(edgeNodeOneXZ, edgeNodeTwoXZ, edgeNodeThreeXZ, edgeNodeFourXZ)
if (periodicAxisId(2).and.periodicAxisId(3)) deallocate(edgeNodeOneYZ, edgeNodeTwoYZ, edgeNodeThreeYZ, edgeNodeFourYZ)
if (periodicity == 3) deallocate(cornerNodeOneYY, cornerNodeTwoYY, cornerNodeThreeYY, cornerNodeFourYY, &
                                 cornerNodeOneXX, cornerNodeTwoXX, cornerNodeThreeXX, cornerNodeFourXX)
deallocate(numElementsOfNode)
deallocate(globalNodeIdTypeDomain)
deallocate(dsEdwMatrix, xsEdwMatrix, coeffEdwMatrix)
deallocate(qqMatrix, qqMatrixFinal, qqMatrixInterpGrafted)
deallocate(phiMatrix, phiTotal)
if (matrixExist.eq.1) then
  deallocate(qqMatrixInterp, dsConvMatrix, xsConvMatrix, coeffConvMatrix)
endif
if (graftedExist.eq.1) then
  deallocate(dsEdwGrafted, dsConvGrafted, xsEdwGrafted, xsConvGrafted, coeffEdwGrafted, coeffConvGrafted)
  deallocate(qqGrafted, qqGraftedFinal, qqGraftedInterp)
  deallocate(graftPointId, deltaNumerical, graftPointValue)
  deallocate(phiGrafted, phiGraftedIndiv)
  if ((exportPhiIndividual.eq.1).AND.(exportAllGraftedChains.eq.0)) deallocate(graftingPointIndexToExport)
endif
deallocate(wwField, wwFieldNew, wwFieldMixed, uuField)
deallocate(nodeVolume)
deallocate(planarCellId, distanceFromFace, planarCellVolume)
deallocate(sphericalCellId, distanceFromNanop, sphericalCellVolume)
deallocate(nodePairId)
deallocate(elementOfNode)
deallocate(nodeBelongsToDirichletFace, nodeBelongsToFaceId)
deallocate(rdiag1)
deallocate(F_m%row, F_m%col, F_m%g, F_m%rh, F_m%c, F_m%k, F_m%w, F_m%isZero)

#ifdef USE_MPI
! Root will send a stop signal to the slave processes
if (root) then
  flag_continue = .False.
  CALL MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
end if

1000 CALL MPI_FINALIZE(ierr)
#endif
!------------------------------------------------------------------------------------------------------------------!
end program RuSseL
