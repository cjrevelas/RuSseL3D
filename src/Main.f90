!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

program RuSseL
!----------------------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod
use arrays_mod
use eos_mod
use constants_mod
use error_handing_mod
use write_helper_mod
use geometry_mod
use iofiles_mod
use delta_mod
use hist_mod
use kcw_mod, only: rdiag1, F_m
use flags_mod, only: contour_symm, contour_asymm, contour_hybrid
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
integer :: ii, kk, iter, ToolsSystemTime, tInit, tFinal, graftNodeId

logical :: convergence = .false., computeDelta = .false.

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
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

if (my_id==0) then
  root = .true.
else
  root = .false.
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
  call MPI_BCAST(mumpsMatrixType, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
  do while (.true.)
    call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (flag_continue) then
      call SolverMumps(mumpsMatrixType)
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
call ParserInput()
call ParserMesh(elemcon)
call InitArrays()

do ii = 1, numNodes
  call ComputeNodeVolume(nodeVolume(ii), ii)
enddo

allocate(dphi2_dr2(numNodes))
dphi2_dr2=0.0d0

call InitDelta()
call ToolsHistogram(binThickness, nodeVolume)

#ifdef USE_MPI
call MPI_BCAST(mumpsMatrixType, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif

call InitChainContour(contourMatrix, lengthMatrixMax, critContourMatrix, numEdwPointsMatrix, stepEdwAveMatrix, ds_mx_ed, xs_mx_ed, coeff_mx_ed)

! In absence of mx chains, skip discretization for mx convolution
if (matrixExist.eq.1) then
  if (contourMatrix.ne.contour_uniform) then
    call InitChainContour(contour_symm, lengthMatrix, critContourMatrix, numConvolPointsMatrix, stepConvolAveMatrix, ds_mx_conv, xs_mx_conv, coeff_mx_conv)
  else
    call InitChainContour(contourMatrix, lengthMatrix, critContourMatrix, numConvolPointsMatrix, stepConvolAveMatrix, ds_mx_conv, xs_mx_conv, coeff_mx_conv)
  endif
endif

if (graftedExist.eq.1) then
  call InitChainContour(contourGrafted, lengthGrafted, critContourGrafted, numEdwPointsGrafted, stepEdwAveGrafted, ds_gr_ed, xs_gr_ed, coeff_gr_ed)
  if (contourGrafted.ne.contour_uniform) then
    call InitChainContour(contour_symm, lengthGrafted, critContourGrafted, numConvolPointsGrafted, stepConvolAveGrafted, ds_gr_conv, xs_gr_conv, coeff_gr_conv)
  else
    call InitChainContour(contourGrafted, lengthGrafted, critContourGrafted, numConvolPointsGrafted, stepConvolAveGrafted, ds_gr_conv, xs_gr_conv, coeff_gr_conv)
  endif
endif

call InitField(uuField, wwField)

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

  call FemMatrixAssemble(rg2OfMatrixMonomer, wwField)

  do ii = 1, numNodes
    qqMatrix(1,ii)      = 1.0d0
    qqMatrixFinal(1,ii) = 1.0d0
  enddo

  ! We always need to solve for mx chains to perform convolution
  call SolverEdwards(ds_mx_ed, numEdwPointsMatrix, qqMatrix, qqMatrixFinal, nodeBelongsToDirichletFace, elemcon)

  if (graftedExist.eq.1) then
    do ii = 1, numNodes
      call interp_linear(1, numEdwPointsMatrix+1, xs_mx_ed, qqMatrixFinal(:,ii), numConvolPointsGrafted+1, xs_gr_conv, qqMatrixInterpGrafted(:,ii))
    enddo

    ! Recompute the delta functions if necessary
    if (getICfromDelta.eq.1) then
      computeDelta = ((iter==0) .OR. ((freeEnergyError <= freeEnergyTolForDelta) .AND. (numGraftedChainsError > numGraftedChainsTol)))

      if (computeDelta) then
        call ComputeDeltaNumerical(elemcon, qqMatrixInterpGrafted, wwFieldMixed, deltaNumerical)

        call ExportDelta(qqMatrixInterpGrafted, deltaNumerical, graftPointValue)
      endif

      do ii = 1, targetNumGraftedChains
        graftNodeId = graftPointId(ii)

        graftPointValue(ii) = deltaNumerical(ii) * lengthGrafted * 1.0d0 / &
                              (qqMatrixInterpGrafted(numConvolPointsGrafted+1,graftNodeId) * (molarBulkDensity * n_avog))
      enddo
    endif

    call FemMatrixAssemble(rg2OfGraftedMonomer, wwField)

    qqGrafted      = 0.0d0
    qqGraftedFinal = 0.0d0

    do ii = 1, targetNumGraftedChains
      graftNodeId = graftPointId(ii)

      qqGrafted(1,graftNodeId)      = graftPointValue(ii)
      qqGraftedFinal(1,graftNodeId) = graftPointValue(ii)
    enddo

    call SolverEdwards(ds_gr_ed, numEdwPointsGrafted, qqGrafted, qqGraftedFinal, nodeBelongsToDirichletFace, elemcon)
  endif

  ! In absence of mx chains, skip mx convolution
  if (matrixExist.eq.1) then
    do ii = 1, numNodes
      call interp_linear(1, numEdwPointsMatrix+1, xs_mx_ed, qqMatrixFinal(:,ii), numConvolPointsMatrix+1, xs_mx_conv, qqMatrixInterp(:,ii))
    enddo

    call ContourConvolution(lengthMatrix, numConvolPointsMatrix, coeff_mx_conv, qqMatrixInterp, qqMatrixInterp, phiMatrix)
  endif

  if (graftedExist.eq.1) then
    do ii = 1, numNodes
      call interp_linear(1, numEdwPointsGrafted+1, xs_gr_ed, qqGraftedFinal(:,ii), numConvolPointsGrafted+1, xs_gr_conv, qqGraftedInterp(:,ii))
    enddo

    call ContourConvolution(lengthGrafted, numConvolPointsGrafted, coeff_gr_conv, qqGraftedInterp, qqMatrixInterpGrafted, phiGrafted)
  endif

  phiTotal = 0.0d0
  do kk = 1, numNodes
    if (matrixExist.eq.1)  phiTotal(kk) = phiTotal(kk) + phiMatrix(kk)
    if (graftedExist.eq.1) phiTotal(kk) = phiTotal(kk) + phiGrafted(kk)
  enddo

  if (matrixExist.eq.1)  call ComputePartitionMatrix(numNodes, numConvolPointsMatrix, qqMatrixInterp, partitionMatrixChains)
  if (matrixExist.eq.1)  call ComputeNumberOfChains(lengthMatrix, phiMatrix, numMatrixChains)
  if (graftedExist.eq.1) call ComputeNumberOfChains(lengthGrafted, phiGrafted, numGraftedChains)

  ! Compare this to russel1d
  do kk = 1, numNodes
    wwFieldNew(kk) = (eos_df_drho(phiTotal(kk)) - eos_df_drho(1.0d0)) / (boltz_const_Joule_K*temperature) - &
                   & sgtParam * (segmentBulkDensity * dphi2_dr2(kk)) / (boltz_const_Joule_K * temperature) + uuField(kk)
  enddo

  wwFieldError    = 0.0d0
  wwFieldStdError = 0.0d0
  wwFieldMaximum  = 0.0d0

  do kk = 1, numNodes
    wwFieldError    = MAX(wwFieldError,DABS(wwFieldNew(kk) - wwField(kk)))
    wwFieldStdError = wwFieldStdError + (wwFieldNew(kk) - wwField(kk))**2.0d0
    wwFieldMaximum  = MAX(wwFieldMaximum, wwFieldNew(kk))
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

  do kk = 1, numNodes
    wwFieldMixed(kk) = (1.0d0 - frac) * wwField(kk) + frac * wwFieldNew(kk)
  enddo

  convergence = (wwFieldError<=fieldTol).OR.(freeEnergyError<=freeEnergyTol)

  call ExportFieldBinary(wwFieldMixed, numNodes, 0)

  if (export(exportFieldBin, iter, convergence)) call ExportFieldBinary(wwFieldMixed, numNodes, iter)

  if ((MOD(iter,1).eq.0).OR.convergence) call ExportEnergies(qqMatrixInterpGrafted, qqGraftedInterp, phiTotal, wwFieldNew, uuField, partitionMatrixChains, targetNumGraftedChains, graftPointId, freeEnergy)

  call ExportComputes(iter, convergence, elemcon)

  call ExportVtuProfiles(phiMatrix, phiGrafted, wwFieldMixed)

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

write(iow,'(3X,A40,E16.9)')adjl("Free energy (mJ/m2):",40),                 freeEnergy
write(6  ,'(3X,A40,E16.9)')adjl("Free energy (mJ/m2):",40),                 freeEnergy
write(iow,'(3X,A40,E16.9)')adjl("Interface area (A2):",40),                 interfaceArea()
write(6  ,'(3X,A40,E16.9)')adjl("Interface area (A2):",40),                 interfaceArea()
write(iow,'(3X,A40,E16.9)')adjl("Partition function of matrix chains:",40), partitionMatrixChains
write(6  ,'(3X,A40,E16.9)')adjl("Partition function of matrix chains:",40), partitionMatrixChains
write(iow,'(3X,A40,E16.9)')adjl("Grafting density (A^-2):",40),             numGraftedChains/interfaceArea()
write(6  ,'(3X,A40,E16.9)')adjl("Grafting density (A^-2):",40),             numGraftedChains/interfaceArea()
write(iow,'(3X,A40,E16.4)')adjl("Number of grafted chains:",40),            numGraftedChains
write(6  ,'(3X,A40,E16.4)')adjl("Number of grafted chains:",40),            numGraftedChains
write(iow,'(3X,A40,E16.4)')adjl("Number of matrix chains:",40),             numMatrixChains
write(6  ,'(3X,A40,E16.4)')adjl("Number of matrix chains:",40),             numMatrixChains

tFinal = ToolsSystemTime()
write(6,'(3X,A40,I16)')adjl('Run duration:',40), tFinal - tInit

! Deallocate all remaining dynamic memory
deallocate(nodeCoord)
deallocate(dphi2_dr2, d2phi_dr2)
if (numDirichletFaces > 0)    deallocate(dirichletFaceId, dirichletFaceValue, plateAlpha, plateSigma)
if (numNanopFaces > 0) deallocate(nanopFaceId, nanopFaceValue, nanopAlpha, nanopSigma, nanopRadiusEff, nanopCenter)
deallocate(numElementsOfNode)
deallocate(globalNodeIdTypeDomain)
deallocate(ds_mx_ed, xs_mx_ed, coeff_mx_ed)
deallocate(qqMatrix, qqMatrixFinal, qqMatrixInterpGrafted)
deallocate(phiMatrix, phiTotal)
if (matrixExist.eq.1) then
  deallocate(qqMatrixInterp, ds_mx_conv, xs_mx_conv, coeff_mx_conv)
endif
if (graftedExist.eq.1) then
  deallocate(ds_gr_ed, ds_gr_conv, xs_gr_ed, xs_gr_conv, coeff_gr_ed, coeff_gr_conv)
  deallocate(qqGrafted, qqGraftedFinal, qqGraftedInterp)
  deallocate(graftPointId, deltaNumerical, graftPointValue)
  deallocate(phiGrafted, phiGraftedIndiv)
  if ((exportPhiIndividual.eq.1).AND.(exportAllGraftedChains.eq.0)) deallocate(gpIndexToExport)
endif
deallocate(wwField, wwFieldNew, wwFieldMixed, uuField)
deallocate(nodeVolume)
deallocate(planar_cell_of_np, dist_from_face, cell_vol_planar)
deallocate(sph_cell_of_np, dist_from_np, cell_vol_sph)
deallocate(nodePairId)
deallocate(elementOfNode)
deallocate(nodeBelongsToDirichletFace, nodeBelongsToFaceId)
deallocate(rdiag1)
deallocate(F_m%row, F_m%col, F_m%g, F_m%rh, F_m%c, F_m%k, F_m%w, F_m%is_zero)

#ifdef USE_MPI
! Root will send a stop signal to the slave processes
if (root) then
  flag_continue = .false.
  call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
end if

1000 call MPI_FINALIZE(ierr)
#endif
!------------------------------------------------------------------------------------------------------------------!
end program RuSseL
