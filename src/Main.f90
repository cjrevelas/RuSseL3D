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
integer :: ii, kk, iter, ToolsSystemTime, t_init, t_final, gnode_id

logical :: convergence = .false., calc_delta = .false.

real(8), allocatable, dimension(:) :: dphi2_dr2

type(fhash_type__ints_double) :: elemcon

real(8) :: partitionMatrixChains = 0.0d0, numMatrixChains = 0.0d0
real(8) :: numGraftedChains = 0.0d0, numGraftedChainsError = 1.0d2
real(8) :: fieldMaximum = 0.0d0, fieldStdError = 0.0d0, fieldError = 2.0d5
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
  call ComputeNodeVolume(volnp(ii), ii)
enddo

allocate(dphi2_dr2(numNodes))
dphi2_dr2=0.0d0

call InitDelta()
call ToolsHistogram(binThickness, volnp)

#ifdef USE_MPI
call MPI_BCAST(mumpsMatrixType, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif

if (matrixExist.eq.1) then
  call InitChainContour(contourMatrix, lengthMatrixMax, critContourMatrix, numEdwPointsMatrix, stepEdwAveMatrix, ds_mx_ed, xs_mx_ed, coeff_mx_ed)
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

call InitField(Ufield, ww)

ww_mix = ww
!**************************************************************************************************************!
!                                        LOOPS FOR FIELD CONVERGENCE                                           !
!**************************************************************************************************************!
write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-----------------------------------SIMULATION STARTS-----------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------SIMULATION STARTS-----------------------------------',85)

write(iow,'(A10,1X,8(A19,1X),A16)') "iter", "fraction", "energy", "energy_error", "gr_chains", "gr_chains_error", "field_error", "field_std_error", "field_max"
write(6  ,'(A4,1X,8(A14,1X),A12)')  "iter", "fraction", "energy", "energy_error", "gr_chains", "gr_chains_error", "field_error", "field_std_error", "field_max"

t_init = ToolsSystemTime()

do iter = initialIterationId, iterations-1
  write(iow,'(I10,1X,8(E19.9E3,1X))') iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, fieldError, fieldStdError, fieldMaximum
  write(6  ,'(I4 ,1X,8(E14.4E3,1X))') iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, fieldError, fieldStdError, fieldMaximum

  close(iow)
  open(unit=iow, file = IO_logFile, position = 'append')

  ww = ww_mix

  call FemMatrixAssemble(rg2OfMatrixMonomer, ww)

  do ii = 1, numNodes
    qmx(1,ii)       = 1.0d0
    qmx_final(1,ii) = 1.0d0
  enddo

  call SolverEdwards(ds_mx_ed, numEdwPointsMatrix, mumpsMatrixType, qmx, qmx_final, nodeBelongsToDirichletFace, elemcon)

  if (graftedExist.eq.1) then
    do ii = 1, numNodes
      call interp_linear(1, numEdwPointsMatrix+1, xs_mx_ed, qmx_final(:,ii), numConvolPointsGrafted+1, xs_gr_conv, qmx_interp_mg(:,ii))
    enddo

    ! Recompute the delta functions if necessary
    if (getICfromDelta.eq.1) then
      calc_delta = ((iter==0) .OR. ((freeEnergyError <= freeEnergyTolForDelta) .AND. (numGraftedChainsError > numGraftedChainsTol)))

      if (calc_delta) then
        call ComputeDeltaNumerical(numNodes, elemcon, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, ww_mix, &
                                   targetNumGraftedChains, graftPointId, deltaNumerical, volnp)
        call ExportDelta(numNodes, qmx_interp_mg, numConvolPointsGrafted, targetNumGraftedChains, graftPointId, deltaNumerical, graftPointValue, volnp)
      endif

      do ii = 1, targetNumGraftedChains
        gnode_id = graftPointId(ii)
        graftPointValue(ii) = deltaNumerical(ii) * lengthGrafted * 1.0d0 / (qmx_interp_mg(numConvolPointsGrafted+1,gnode_id) * (molarBulkDensity * n_avog))
      enddo
    endif

    call FemMatrixAssemble(rg2OfGraftedMonomer, ww)

    qgr       = 0.0d0
    qgr_final = 0.0d0

    do ii = 1, targetNumGraftedChains
      gnode_id = graftPointId(ii)

      qgr(1,gnode_id)       = graftPointValue(ii)
      qgr_final(1,gnode_id) = graftPointValue(ii)
    enddo

    call SolverEdwards(ds_gr_ed, numEdwPointsGrafted, mumpsMatrixType, qgr, qgr_final, nodeBelongsToDirichletFace, elemcon)
  endif

  if (matrixExist.eq.1) then
    do ii = 1, numNodes
      call interp_linear(1, numEdwPointsMatrix+1, xs_mx_ed, qmx_final(:,ii), numConvolPointsMatrix+1, xs_mx_conv, qmx_interp_mm(:,ii))
    enddo

    call ContourConvolution(numNodes, lengthMatrix, numConvolPointsMatrix, coeff_mx_conv, qmx_interp_mm, qmx_interp_mm, phi_mx)
  endif

  if (graftedExist.eq.1) then
    do ii = 1, numNodes
      call interp_linear(1, numEdwPointsGrafted+1, xs_gr_ed, qgr_final(:,ii), numConvolPointsGrafted+1, xs_gr_conv, qgr_interp(:,ii))
    enddo

    call ContourConvolution(numNodes, lengthGrafted, numConvolPointsGrafted, coeff_gr_conv, qgr_interp, qmx_interp_mg, phi_gr)
  endif

  phi_total = 0.0d0
  do kk = 1, numNodes
    if (matrixExist.eq.1) phi_total(kk) = phi_total(kk) + phi_mx(kk)
    if (graftedExist.eq.1) phi_total(kk) = phi_total(kk) + phi_gr(kk)
  enddo

  if (matrixExist.eq.1)  call ComputePartitionMatrix(numNodes, numConvolPointsMatrix, qmx_interp_mm, partitionMatrixChains)
  if (matrixExist.eq.1)  call ComputeNumberOfChains(numNodes, lengthMatrix, molarBulkDensity, phi_mx, numMatrixChains)
  if (graftedExist.eq.1) call ComputeNumberOfChains(numNodes, lengthGrafted, molarBulkDensity, phi_gr, numGraftedChains)

  do kk = 1, numNodes
    ww_new(kk) = (eos_df_drho(phi_total(kk)) - eos_df_drho(1.0d0)) / (boltz_const_Joule_K*temperature) - &
                 & sgtParam * (segmentBulkDensity * dphi2_dr2(kk)) / (boltz_const_Joule_K * temperature) + Ufield(kk)
  enddo

  fieldError    = 0.0d0
  fieldStdError = 0.0d0
  fieldMaximum  = 0.0d0

  do kk = 1, numNodes
    fieldError    = MAX(fieldError,DABS(ww_new(kk) - ww(kk)))
    fieldStdError = fieldStdError + (ww_new(kk) - ww(kk))**2.0d0
    fieldMaximum  = MAX(fieldMaximum, ww_new(kk))
  enddo

  fieldStdError = SQRT(fieldStdError / FLOAT((numNodes - 1)))
  fieldMaximum  = fieldMaximum  * lengthMatrix
  fieldError    = fieldError    * lengthMatrix
  fieldStdError = fieldStdError * lengthMatrix

  freeEnergyError    = ABS(freeEnergy - freeEnergyPrevious)
  freeEnergyPrevious = freeEnergy

  if (graftedExist.eq.1) then
    numGraftedChainsError = ABS(numGraftedChains-DBLE(targetNumGraftedChains)) / DBLE(targetNumGraftedChains)
  else
    numGraftedChainsError = 0.0d0
  endif

  do kk = 1, numNodes
    ww_mix(kk) = (1.0d0 - frac) * ww(kk) + frac * ww_new(kk)
  enddo

  convergence = (fieldError<=fieldTol).OR.(freeEnergyError<=freeEnergyTol)

  call ExportFieldBinary(ww_mix, numNodes, 0)

  if (export(exportFieldBin, iter, convergence)) call ExportFieldBinary(ww_mix, numNodes, iter)

  if ((MOD(iter,1).eq.0).OR.convergence) call ExportEnergies(qmx_interp_mg, qgr_interp, phi_total, ww_new, Ufield, partitionMatrixChains, targetNumGraftedChains, graftPointId, freeEnergy)

  call ExportComputes(iter, convergence, elemcon)

  call ExportVtu(phi_mx)

  if (convergence) exit
enddo
!**************************************************************************************************************!
!                                             EXPORT SIMULATION RESULTS                                        !
!**************************************************************************************************************!
write(iow,'(I10,1X,8(E19.9E3,1X))')  iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, fieldError, fieldStdError, fieldMaximum
write(6  ,'(I4 ,1X,8(E14.4E3,1X))')  iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, fieldError, fieldStdError, fieldMaximum

write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-----------------------------------SUMMARIZED RESULTS-----------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------SUMMARIZED RESULTS-----------------------------------',85)

if (fieldError.lt.fieldTol) then
  write(iow,'("Field convergence of max error",F16.9)') fieldError
  write(6  ,'("Field convergence of max error",F16.9)') fieldError
endif

if (freeEnergyError.lt.freeEnergyTol) then
  write(iow,'("Energy convergence of max error",F16.9)') freeEnergyError
  write(6  ,'("Energy convergence of max error",F16.9)') freeEnergyError
endif

write(iow,'(3X,A40,E16.9)')adjl("Free energy (mJ/m2):",40),                 freeEnergy
write(6  ,'(3X,A40,E16.9)')adjl("Free energy (mJ/m2):",40),                 freeEnergy
write(iow,'(3X,A40,E16.9)')adjl("Interface area (A2):",40),                 interf_area()
write(6  ,'(3X,A40,E16.9)')adjl("Interface area (A2):",40),                 interf_area()
write(iow,'(3X,A40,E16.9)')adjl("Partition function of matrix chains:",40), partitionMatrixChains
write(6  ,'(3X,A40,E16.9)')adjl("Partition function of matrix chains:",40), partitionMatrixChains
write(iow,'(3X,A40,E16.9)')adjl("Grafting density (A^-2):",40),             numGraftedChains/interf_area()
write(6  ,'(3X,A40,E16.9)')adjl("Grafting density (A^-2):",40),             numGraftedChains/interf_area()
write(iow,'(3X,A40,E16.4)')adjl("Number of grafted chains:",40),            numGraftedChains
write(6  ,'(3X,A40,E16.4)')adjl("Number of grafted chains:",40),            numGraftedChains
write(iow,'(3X,A40,E16.4)')adjl("Number of matrix chains:",40),             numMatrixChains
write(6  ,'(3X,A40,E16.4)')adjl("Number of matrix chains:",40),             numMatrixChains

t_final = ToolsSystemTime()
write(6,'(3X,A40,I16)')adjl('Run duration:',40), t_final - t_init

! Deallocate all remaining dynamic memory
deallocate(nodeCoord)
deallocate(dphi2_dr2, d2phi_dr2)
if (numDirichletFaces > 0)    deallocate(dirichletFaceId, dirichletFaceValue, plateAlpha, plateSigma)
if (numNanopFaces > 0) deallocate(nanopFaceId, nanopFaceValue, nanopAlpha, nanopSigma, nanopRadiusEff, nanopCenter)
deallocate(numElementsOfNode)
deallocate(globalNodeIdTypeDomain)
deallocate(ds_mx_ed, xs_mx_ed, coeff_mx_ed)
deallocate(qmx, qmx_final, qmx_interp_mg)
deallocate(phi_mx, phi_total)
if (matrixExist.eq.1) then
  deallocate(qmx_interp_mm, ds_mx_conv, xs_mx_conv, coeff_mx_conv)
endif
if (graftedExist.eq.1) then
  deallocate(ds_gr_ed, ds_gr_conv, xs_gr_ed, xs_gr_conv, coeff_gr_ed, coeff_gr_conv)
  deallocate(qgr, qgr_final, qgr_interp)
  deallocate(graftPointId, deltaNumerical, graftPointValue)
  deallocate(phi_gr, phi_gr_indiv)
endif
deallocate(ww, ww_new, ww_mix, Ufield)
deallocate(volnp)
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
