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
integer :: ii, kk, iter, tools_sys_time, t_init, t_final, gnode_id

logical :: convergence = .false., calc_delta = .false.

real(8), allocatable, dimension(:) :: dphi2_dr2

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
      call solver_mumps(mumpsMatrixType)
    else
      exit
    endif
  end do
  goto 1000
endif
#endif

iow = 10
open(unit=iow, file = logfile)

! Initialize the error log
ioe = 11
open(unit=ioe, file = errorfile, status='replace')
close(ioe)
!**************************************************************************************************************!
!                                             INITIALIZATION SECTION                                           !
!**************************************************************************************************************!
call parser_input
call parser_mesh
call init_arrays


do ii = 1, numnp
  call compute_node_volume(volnp(ii), ii)
enddo

allocate(dphi2_dr2(numnp))
dphi2_dr2=0.0d0

call init_delta
call tools_histogram(binThickness, volnp)

#ifdef USE_MPI
call MPI_BCAST(mumpsMatrixType, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif

if (mx_exist.eq.1) then
  call init_chain_contour(contour_discr_mx, lengthMatrixMax, xs_crit_mx, ns_mx_ed, ds_ave_mx_ed, ds_mx_ed, xs_mx_ed, coeff_mx_ed)
  if (contour_discr_mx.ne.contour_uniform) then
    call init_chain_contour(contour_symm, lengthMatrix, xs_crit_mx, ns_mx_conv, ds_ave_mx_conv, ds_mx_conv, xs_mx_conv, coeff_mx_conv)
  else
    call init_chain_contour(contour_discr_mx, lengthMatrix, xs_crit_mx, ns_mx_conv, ds_ave_mx_conv, ds_mx_conv, xs_mx_conv, coeff_mx_conv)
  endif
endif

if (gr_exist.eq.1) then
  call init_chain_contour(contour_discr_gr, lengthGrafted, xs_crit_gr, ns_gr_ed, ds_ave_gr_ed, ds_gr_ed, xs_gr_ed, coeff_gr_ed)
  if (contour_discr_gr.ne.contour_uniform) then
    call init_chain_contour(contour_symm, lengthGrafted, xs_crit_gr, ns_gr_conv, ds_ave_gr_conv, ds_gr_conv, xs_gr_conv, coeff_gr_conv)
  else
    call init_chain_contour(contour_discr_gr, lengthGrafted, xs_crit_gr, ns_gr_conv, ds_ave_gr_conv, ds_gr_conv, xs_gr_conv, coeff_gr_conv)
  endif
endif

call init_field(Ufield, ww)

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

t_init = tools_sys_time()

do iter = init_iter, iterations-1
    write(iow,'(I10,1X,8(E19.9E3,1X))') iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, fieldError, fieldStdError, fieldMaximum
    write(6  ,'(I4 ,1X,8(E14.4E3,1X))') iter, frac, freeEnergy, freeEnergyError, numGraftedChains, numGraftedChainsError, fieldError, fieldStdError, fieldMaximum

    close(iow)
    open(unit=iow, file = logfile, position = 'append')

    ww = ww_mix

    call fem_matrix_assemble(rg2OfMatrixMonomer, ww)

    do ii = 1, numnp
        qmx(1,ii)       = 1.0d0
        qmx_final(1,ii) = 1.0d0
    enddo

    call solver_edwards(ds_mx_ed, ns_mx_ed, mumpsMatrixType, qmx, qmx_final, nodeBelongsToDirichletFace)

    if (gr_exist.eq.1) then
        do ii = 1, numnp
            call interp_linear(1, ns_mx_ed+1, xs_mx_ed, qmx_final(:,ii), ns_gr_conv+1, xs_gr_conv, qmx_interp_mg(:,ii))
        enddo

        ! Recompute the delta functions if necessary
        if (grafted_ic_from_delta.eq.1) then
            calc_delta = ((iter==0) .OR. ((freeEnergyError <= freeEnergyTolForDelta) .AND. (numGraftedChainsError > numGraftedChainsTol)))

            if (calc_delta) then
                call compute_delta_numer(numnp, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, ww_mix, targetNumGraftedChains, gpid, delta_numer, volnp)
                call export_delta(numnp, qmx_interp_mg, ns_gr_conv, targetNumGraftedChains, gpid, delta_numer, gp_init_value, volnp)
            endif

            do ii = 1, targetNumGraftedChains
                gnode_id = gpid(ii)
                gp_init_value(ii) = delta_numer(ii) * lengthGrafted * 1.0d0 / (qmx_interp_mg(ns_gr_conv+1,gnode_id) * (molarBulkDensity * n_avog))
            enddo
        endif

        call fem_matrix_assemble(rg2OfGraftedMonomer, ww)

        qgr       = 0.0d0
        qgr_final = 0.0d0

        do ii = 1, targetNumGraftedChains
            gnode_id = gpid(ii)

            qgr(1,gnode_id)       = gp_init_value(ii)
            qgr_final(1,gnode_id) = gp_init_value(ii)
        enddo

        call solver_edwards(ds_gr_ed, ns_gr_ed, mumpsMatrixType, qgr, qgr_final, nodeBelongsToDirichletFace)
    endif

    if (mx_exist.eq.1) then
        do ii = 1, numnp
            call interp_linear(1, ns_mx_ed+1, xs_mx_ed, qmx_final(:,ii), ns_mx_conv+1, xs_mx_conv, qmx_interp_mm(:,ii))
        enddo

        call contour_convolution(numnp, lengthMatrix, ns_mx_conv, coeff_mx_conv, qmx_interp_mm, qmx_interp_mm, phi_mx)
    endif

    if (gr_exist.eq.1) then
        do ii = 1, numnp
            call interp_linear(1, ns_gr_ed+1, xs_gr_ed, qgr_final(:,ii), ns_gr_conv+1, xs_gr_conv, qgr_interp(:,ii))
        enddo

        call contour_convolution(numnp, lengthGrafted, ns_gr_conv, coeff_gr_conv, qgr_interp, qmx_interp_mg, phi_gr)
    endif

    phi_total = 0.0d0
    do kk = 1, numnp
        if (mx_exist.eq.1) phi_total(kk) = phi_total(kk) + phi_mx(kk)
        if (gr_exist.eq.1) phi_total(kk) = phi_total(kk) + phi_gr(kk)
    enddo

    if (mx_exist.eq.1) call compute_part_func_mx(numnp, ns_mx_conv, qmx_interp_mm, partitionMatrixChains)
    if (mx_exist.eq.1) call compute_number_of_chains(numnp, lengthMatrix, molarBulkDensity, phi_mx, numMatrixChains)
    if (gr_exist.eq.1) call compute_number_of_chains(numnp, lengthGrafted, molarBulkDensity, phi_gr, numGraftedChains)

    do kk = 1, numnp
        ww_new(kk) = (eos_df_drho(phi_total(kk)) - eos_df_drho(1.0d0)) / (boltz_const_Joule_K*temperature) - &
                   & k_gr * (segmentBulkDensity * dphi2_dr2(kk)) / (boltz_const_Joule_K * temperature) + Ufield(kk)
    enddo

    fieldError    = 0.0d0
    fieldStdError = 0.0d0
    fieldMaximum  = 0.0d0

    do kk = 1, numnp
        fieldError    = MAX(fieldError,DABS(ww_new(kk) - ww(kk)))
        fieldStdError = fieldStdError + (ww_new(kk) - ww(kk))**2.0d0
        fieldMaximum  = MAX(fieldMaximum, ww_new(kk))
    enddo

    fieldStdError = SQRT(fieldStdError / FLOAT((numnp - 1)))
    fieldMaximum       = fieldMaximum       * lengthMatrix
    fieldError  = fieldError  * lengthMatrix
    fieldStdError = fieldStdError * lengthMatrix

    freeEnergyError = ABS(freeEnergy - freeEnergyPrevious)
    freeEnergyPrevious  = freeEnergy

    if (gr_exist.eq.1) then
        numGraftedChainsError = ABS(numGraftedChains-DBLE(targetNumGraftedChains)) / DBLE(targetNumGraftedChains)
    else
        numGraftedChainsError = 0.0d0
    endif

    do kk = 1, numnp
        ww_mix(kk) = (1.0d0 - frac) * ww(kk) + frac * ww_new(kk)
    enddo

    convergence = (fieldError<=fieldTol).OR.(freeEnergyError<=freeEnergyTol)

    call export_field_bin(ww_mix, numnp, 0)

    if (export(exportFieldBinary, iter, convergence)) call export_field_bin(ww_mix, numnp, iter)

    if ((MOD(iter,1).eq.0).OR.convergence) call export_energies(qmx_interp_mg, qgr_interp, phi_total, ww_new, Ufield, partitionMatrixChains, targetNumGraftedChains, gpid, freeEnergy)

    call export_computes(iter, convergence)

    call export_vtu(phi_mx)

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

t_final = tools_sys_time()
write(6,'(3X,A40,I16)')adjl('Run duration:',40), t_final - t_init

#ifdef USE_MPI
! Root will send a stop signal to the slave processes
if (root) then
    flag_continue = .false.
    call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
end if

1000 call MPI_FINALIZE(ierr)
#endif

! Deallocate all remaining dynamic memory
deallocate(xc)
deallocate(dphi2_dr2, d2phi_dr2)
if (numDirichletFaces > 0)    deallocate(dirichletFaceId, A_plate, sigma_plate)
if (numNanoparticleFaces > 0) deallocate(nanoparticleFaceId, A_np, sigma_np, radius_np_eff, center_np)
deallocate(num_of_elems_of_node)
deallocate(global_node_id_type_domain)
deallocate(ds_mx_ed, xs_mx_ed, coeff_mx_ed)
deallocate(qmx, qmx_final, qmx_interp_mg)
deallocate(phi_mx, phi_total)
if (mx_exist.eq.1) then
    deallocate(qmx_interp_mm, ds_mx_conv, xs_mx_conv, coeff_mx_conv)
endif
if (gr_exist.eq.1) then
    deallocate(ds_gr_ed, ds_gr_conv, xs_gr_ed, xs_gr_conv, coeff_gr_ed, coeff_gr_conv)
    deallocate(qgr, qgr_final, qgr_interp)
    deallocate(gpid, delta_numer, gp_init_value)
    deallocate(phi_gr, phi_gr_indiv)
endif
deallocate(ww, ww_new, ww_mix, Ufield)
deallocate(volnp)
deallocate(planar_cell_of_np, dist_from_face, cell_vol_planar)
deallocate(sph_cell_of_np, dist_from_np, cell_vol_sph)
deallocate(node_pair_id)
deallocate(el_node)
deallocate(nodeBelongsToDirichletFace)
deallocate(rdiag1)
deallocate(F_m%row, F_m%col, F_m%g, F_m%rh, F_m%c, F_m%k, F_m%w, F_m%is_zero)
!------------------------------------------------------------------------------------------------------------------!
end program RuSseL
