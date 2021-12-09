!----------------------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!
!                 THIS IS THE FINITE ELEMENT CODE WHICH RUNS A SELF-CONSISTENT FIELD SIMULATION IN THREE DIMENSIONS.               !
!          IT APPLIES AN ANDERSON MIXING RULE WITH RESPECT TO THE FIELD IN ORDER TO CONVERGE TO A SELF-CONSISTENT SOLUTION.        !
!      IT OFFERS THE ABILITY TO WORK WITH BOTH MATRIX AND GRAFTED CHAINS IN THE PRESENCE OF PLANAR SOLID SURFACES OR NANOPARTICLES.!
!                                                      DATE: 01/07/2019                                                            !
!                         AUTHORS: CONSTANTINOS J. REVELAS, ARISTOTELIS P. SGOUROS, APOSTOLIS T. LAKKAS                            !
!----------------------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!
program FEM_3D
!----------------------------------------------------------------------------------------------------------------------------------!
use parser_vars
use arrays
use eos
use constants
use error_handing
use write_helper
use geometry
use iofiles
use delta
use flags, only: symm, asymm, hybrid
#ifdef USE_MPI
use mpistuff
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
#ifdef USE_MPI
include "mpif.h"
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
integer :: ii, kk, iter, get_sys_time, t_init, t_final, gnode_id!, node_id

logical :: convergence = .false.

!real(8) :: compute_gradient, interp_fem
!real(8) :: phi_plus_dx = 0.d0, phi_plus_dy = 0.d0, phi_plus_dz = 0.d0
!real(8) :: phi_minus_dx = 0.d0, phi_minus_dy = 0.d0, phi_minus_dz = 0.d0
!real(8) :: dphi2_dx2 = 0.d0, dphi2_dy2 = 0.d0, dphi2_dz2 = 0.d0
real(8), allocatable, dimension(:) :: dphi2_dr2

real(8) :: part_func = 0.d0, nch_gr = 0.d0
real(8) :: wa_max = 0.d0, wa_std_error = 0.d0, max_error = 200000.d0
real(8) :: free_energy_prev = 1.d10, free_energy = 0.d0, free_energy_error = 1.d10, free_energy_error_tol = 1.d-6
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

!the slaves will enter the mumps subroutine until they receive a stop signal from master proc
if (.not.root) then
    !receive the matrix type from root
    call MPI_BCAST(mumps_matrix_type, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    do while (.true.)
        call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if (flag_continue) then
            call mumps_sub(mumps_matrix_type)
        else
            exit
        endif
    end do
    goto 1000
endif
#endif

iow = 10
open(unit=iow, file = logfile)

!initialize the error log
ioe = 11
open(unit=ioe, file = errorfile, status='replace')
close(ioe)
!**************************************************************************************************************!
!                                             INITIALIZATION SECTION                                           !
!**************************************************************************************************************!
call parser
call init_scf_params
if ((mx_exist.eq.1).and.(discr_mx.eq.hybrid)) then
   call get_contour_step(ds_ave_mx_ed, xs_crit_mx, chainlen_mx_max, ns_mx_ed)
endif
if ((gr_exist.eq.1).and.(discr_gr.eq.hybrid)) then
   call get_contour_step(ds_ave_gr_ed, xs_crit_gr, chainlen_gr, ns_gr_ed)
endif
call import_mesh
call init_arrays

do ii = 1, numnp
   call get_volnp(volnp(ii), ii)
enddo

allocate(dphi2_dr2(numnp))
dphi2_dr2=0.d0

call init_delta
call alloc_hist(0.1d0, volnp)

#ifdef USE_MPI
call MPI_BCAST(mumps_matrix_type, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif

if (mx_exist.eq.1) then
    call init_chain_contour(discr_mx, chainlen_mx_max, ns_mx_ed, ds_ave_mx_ed, ds_mx_ed, xs_mx_ed, coeff_mx_ed)
    call init_chain_contour(symm, chainlen_mx, ns_mx_conv, ds_ave_mx_conv, ds_mx_conv, xs_mx_conv, coeff_mx_conv)
endif

if (gr_exist.eq.1) then
    call init_chain_contour(discr_gr, chainlen_gr, ns_gr_ed, ds_ave_gr_ed, ds_gr_ed, xs_gr_ed, coeff_gr_ed)
    call init_chain_contour(symm, chainlen_gr, ns_gr_conv, ds_ave_gr_conv, ds_gr_conv, xs_gr_conv, coeff_gr_conv)
endif

call init_field(Ufield, wa)

wa_mix = wa
!**************************************************************************************************************!
!                                        LOOPS FOR FIELD CONVERGENCE                                           !
!**************************************************************************************************************!
write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-----------------------------------SIMULATION STARTS-----------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------SIMULATION STARTS-----------------------------------',85)

write(iow,'(A10,1X,6(A19,1X),A16)') "iter", "fraction", "free_energy", "n_gr_chains", "max_error", "std_error", "wa_max"
write(6  ,'(A4,1X,6(A14,1X),A12)')  "iter", "fraction", "free_energy", "n_gr_chains", "max_error", "std_error", "wa_max"

t_init = get_sys_time()

do iter = init_iter, iterations-1
    write(iow,'(I10,1X,6(E19.9E3,1X))') iter, frac, free_energy, nch_gr, max_error, wa_std_error, wa_max
    write(6  ,'(I4 ,1X,6(E14.4E3,1X))') iter, frac, free_energy, nch_gr, max_error, wa_std_error, wa_max

    close(iow)
    open(unit=iow, file = logfile, position = 'append')

    wa = wa_mix

    call matrix_assemble(Rg2_per_mon_mx, wa)

    do ii = 1, numnp
        qmx(1,ii)       = 1.d0
        qmx_final(1,ii) = 1.d0
    enddo

    call edwards(ds_mx_ed, ns_mx_ed, mumps_matrix_type, qmx, qmx_final, node_in_q0_face)

    if (gr_exist.eq.1) then
        do ii = 1, numnp
            call interp_linear(1, ns_mx_ed+1, xs_mx_ed, qmx_final(:,ii), ns_gr_conv+1, xs_gr_conv, qmx_interp_mg(:,ii))
        enddo

        ! recompute the delta functions if they exceed
        if (grafted_ic_from_delta.eq.1) then
            if ( (iter.eq.1) .or. (MOD(iter,calc_delta_every).eq.0 .and. (ABS(nch_gr-DBLE(num_gpoints))/DBLE(num_gpoints))>0.005d0) ) then
                call compute_delta_numer(numnp, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, wa_mix, num_gpoints, &
     &                             gpid, delta_numer, volnp)
                call export_delta(numnp, qmx_interp_mg, ns_gr_conv, num_gpoints, gpid, delta_numer, gp_init_value, volnp)
            endif

            do ii = 1, num_gpoints
                gnode_id = gpid(ii)
                gp_init_value(ii) = delta_numer(ii) * chainlen_gr &
                                  * 1.d0 / (qmx_interp_mg(ns_gr_conv+1,gnode_id) * (rho_mol_bulk * n_avog))
            enddo
        endif

        call matrix_assemble(Rg2_per_mon_gr, wa)

        qgr       = 0.d0
        qgr_final = 0.d0

        do ii = 1, num_gpoints
            gnode_id = gpid(ii)

            qgr(1,gnode_id)       = gp_init_value(ii)
            qgr_final(1,gnode_id) = gp_init_value(ii)
        enddo

        call edwards(ds_gr_ed, ns_gr_ed, mumps_matrix_type, qgr, qgr_final, node_in_q0_face)
    endif

    if (mx_exist.eq.1) then
        do ii = 1, numnp
            call interp_linear(1, ns_mx_ed+1, xs_mx_ed, qmx_final(:,ii), ns_mx_conv+1, xs_mx_conv, qmx_interp_mm(:,ii))
        enddo

        call convolution(numnp, chainlen_mx, ns_mx_conv, coeff_mx_conv, qmx_interp_mm, qmx_interp_mm, phia_mx)
    endif

    if (gr_exist.eq.1) then
        do ii = 1, numnp
            call interp_linear(1, ns_gr_ed+1, xs_gr_ed, qgr_final(:,ii), ns_gr_conv+1, xs_gr_conv, qgr_interp(:,ii))
        enddo

        call convolution(numnp, chainlen_gr, ns_gr_conv, coeff_gr_conv, qgr_interp, qmx_interp_mg, phia_gr)
    endif

    phi_total = 0.d0
    do kk = 1, numnp
        if (mx_exist.eq.1) phi_total(kk) = phi_total(kk) + phia_mx(kk)
        if (gr_exist.eq.1) phi_total(kk) = phi_total(kk) + phia_gr(kk)
    enddo

    if (mx_exist.eq.1) call get_part_func(numnp, ns_mx_conv, qmx_interp_mm, part_func)

    if (gr_exist.eq.1) then
        call get_nchains(numnp, chainlen_gr, rho_mol_bulk, phia_gr, nch_gr)
    endif

    !do kk = 1, numnp
    !    dphi2_dr2(kk) = compute_gradient(kk, phi_total)
        !STOP
    !enddo
    !STOP
    !phi_plus_dx = interp_fem(1, xc(1,1)+dx, xc(2,1), xc(3,1), phi_total(1))
    !phi_plus_dy = interp_fem(1, xc(1,1), xc(2,1)+dy, xc(3,1), phi_total(1))
    !phi_plus_dz = interp_fem(1, xc(1,1), xc(2,1), xc(3,1)+dz, phi_total(1))

    !dphi2_dx2 = (phi_plus_dx -2.D0*phi_total(1) + phi_total(1)) / (dx*1.D-10)**2.D0
    !dphi2_dy2 = (phi_plus_dy -2.D0*phi_total(1) + phi_total(1)) / (dy*1.D-10)**2.D0
    !dphi2_dz2 = (phi_plus_dz -2.D0*phi_total(1) + phi_total(1)) / (dz*1.D-10)**2.D0

    !dphi2_dr2(1) = dphi2_dx2 + dphi2_dy2 + dphi2_dz2

    !phi_minus_dx = interp_fem(numnp, xc(1,numnp)-dx, xc(2,numnp), xc(3,numnp), phi_total(numnp))
    !phi_minus_dy = interp_fem(numnp, xc(1,numnp), xc(2,numnp)-dy, xc(3,numnp), phi_total(numnp))
    !phi_minus_dz = interp_fem(numnp, xc(1,numnp), xc(2,numnp), xc(3,numnp)-dz, phi_total(numnp))

    !dphi2_dx2 = (phi_total(numnp) -2.D0*phi_total(numnp) + phi_minus_dx) / (dx*1.D-10)**2.D0
    !dphi2_dy2 = (phi_total(numnp) -2.D0*phi_total(numnp) + phi_minus_dy) / (dy*1.D-10)**2.D0
    !dphi2_dz2 = (phi_total(numnp) -2.D0*phi_total(numnp) + phi_minus_dz) / (dz*1.D-10)**2.D0

    !dphi2_dr2(numnp) = dphi2_dx2 + dphi2_dy2 + dphi2_dz2

    !do kk = 1, numnp
    !    write(123,*) xc(3,kk), phi_total(kk), dphi2_dr2(kk)
    !enddo
    !node_id = 200
    !dphi2_dr2 = compute_gradient(node_id, phi_total)

    do kk = 1, numnp
        wa_new(kk) = (eos_df_drho(phi_total(kk)) - eos_df_drho(1.d0)) / (boltz_const_Joule_K*Temp) - &
                   & k_gr * (rho_seg_bulk * dphi2_dr2(kk)) / (boltz_const_Joule_K * Temp) + Ufield(kk)
    enddo

    max_error    = 0.d0
    wa_std_error = 0.d0
    wa_max       = 0.d0

    do kk = 1, numnp
        max_error    = MAX(max_error,DABS(wa_new(kk) - wa(kk)))
        wa_std_error = wa_std_error + (wa_new(kk) - wa(kk))**2
        wa_max       = MAX(wa_max, wa_new(kk))
    enddo

    wa_std_error = SQRT(wa_std_error / FLOAT((numnp - 1)))
    wa_max       = wa_max       * chainlen_mx
    max_error    = max_error    * chainlen_mx
    wa_std_error = wa_std_error * chainlen_mx

    free_energy_error = ABS(free_energy - free_energy_prev)
    free_energy_prev  = free_energy

    do kk = 1, numnp
        wa_mix(kk) = (1.d0 - frac) * wa(kk) + frac * wa_new(kk)
    enddo

    if (MOD(iter,output_every).eq.0) call export_field_bin(wa_mix, numnp, iter)

    convergence = (max_error<=max_error_tol).or.(free_energy_error<=free_energy_error_tol)

    if (MOD(iter,1).eq.0.or.convergence) call energies(qmx_interp_mg, phi_total, part_func, num_gpoints, gpid, free_energy)
    if (output_every.gt.0)then
        if (MOD(iter,output_every).eq.0.or.convergence) call export_computes()
    endif
    if (convergence) exit
enddo
!**************************************************************************************************************!
!                                             EXPORT SIMULATION RESULTS                                        !
!**************************************************************************************************************!
call export_field_bin(wa_mix, numnp, iter)

write(iow,'(I10,1X,6(E19.9E3,1X))')  iter, frac, free_energy, nch_gr, max_error, wa_std_error, wa_max
write(6  ,'(I4 ,1X,6(E14.4E3,1X))')  iter, frac, free_energy, nch_gr, max_error, wa_std_error, wa_max

write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-----------------------------------SUMMARIZED RESULTS-----------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------SUMMARIZED RESULTS-----------------------------------',85)

if (max_error.lt.max_error_tol) then
    write(iow,'("Convergence of max error",F16.9)') max_error
    write(6  ,'("Convergence of max error",F16.9)') max_error
endif

write(iow,'(3X,A40,E16.9)')adjl("Free energy (mJ/m2):",40),                 free_energy
write(6  ,'(3X,A40,E16.9)')adjl("Free energy (mJ/m2):",40),                 free_energy
write(iow,'(3X,A40,E16.9)')adjl("Partition function of matrix chains:",40), part_func
write(6  ,'(3X,A40,E16.9)')adjl("Partition function of matrix chains:",40), part_func
write(iow,'(3X,A40,E16.9)')adjl("Grafting density (A^-2):",40),             nch_gr/interf_area
write(6  ,'(3X,A40,E16.9)')adjl("Grafting density (A^-2):",40),             nch_gr/interf_area
write(iow,'(3X,A40,E16.4)')adjl("Number of grafted chains:",40),            nch_gr
write(6  ,'(3X,A40,E16.4)')adjl("Number of grafted chains:",40),            nch_gr

t_final = get_sys_time()
write(6,'(3X,A40,I16)')adjl('Run duration:',40), t_final - t_init

#ifdef USE_MPI
!root will send a stop signal to the slaves
if (root) then
    flag_continue = .false.
    call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
end if

1000 call MPI_FINALIZE(ierr)
#endif
!------------------------------------------------------------------------------------------------------------------!
end program FEM_3D
