program FEM_3D
!-------------------------------------------------------------------!
use parser_vars
use constants
use error_handing
use write_helper
use kcw
use geometry
#ifdef USE_MPI
use mpistuff
#endif
!-------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------!
#ifdef USE_MPI
include 'mpif.h'
#endif
!-------------------------------------------------------------------!
integer :: i1, k1, iter, n_outside
integer :: get_sys_time, t_init, t_final

character(20) :: gp_filename = "gnodes.lammpstrj"
character(20) :: field_in_filename = 'field_in.bin'

real(8)                              :: wa_max = 0.d0, wa_max_abs = 0.d0
real(8)                              :: wa_std_error = 0.d0, max_error = 200000.d0
real(8)                              :: mix_tol = 0.d0, wa_step = 0.d0, wa_ave = 0.d0
real(8)                              :: part_func = 0.d0, nch_gr = 0.d0
real(8)                              :: adh_ten = 0.d0
real(8), allocatable, dimension(:)   :: ds_free, ds_gr
real(8), allocatable, dimension(:)   :: koeff_free, koeff_gr
real(8), allocatable, dimension(:)   :: wa, wa_new, wa_mix, Ufield
real(8), allocatable, dimension(:)   :: phia_fr, phia_gr
real(8), allocatable, dimension(:,:) :: qf, qgr
real(8), allocatable, dimension(:,:) :: qf_final, qgr_final
!-------------------------------------------------------------------!
!*******************************************************************!
!                           MPI SECTION                             !
!*******************************************************************!
#ifdef USE_MPI
call mpi_init(ierr)
call MPI_Comm_size(MPI_COMM_WORLD, n_proc, ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)

if (my_id==0) then
    root = .true.
else
    root = .false.
endif

if (root) then
    write (*,*)
    write (*,'(a,i4,a)') ' MPI run with ', n_proc, ' procs.'
    write (*,*)
end if

flag_continue = .true.

!the slaves will enter the mumps subroutine until they receive a stop
!signal from master proc
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
!-------------------------------------------------------------------!
iow = 10
open(unit=iow, file = 'scft.out.txt')
write(iow,'(''Polymer FEM_3D 1 Sept 19 '')')
write(iow,'(''Polumer Bulk with solid surfaces   ''/    &
          & ''---------------------------------  ''/    &
          & ''SCF theory, Gaussian string model  ''/    &
          & ''Finite Element solution with       ''/    &
          & ''successive substitutions           '')')

!initialize the error log
ioe = 11
open(unit=ioe, file = 'error.out.txt', status='replace')
close(ioe)
!*******************************************************************!
!                       INITIALIZATION SECTION                      !
!*******************************************************************!
!parse the input data from the datafile
call parser

!calculate essential scf parameters
call calc_scf_params

!read the input from the mesh file and generate it
call mesh

!allocate and initialize essential scf arrays
allocate(rdiag1(numnp))
allocate(wa(numnp),wa_mix(numnp),wa_new(numnp),Ufield(numnp))

wa         = 0.d0
wa_mix     = 0.d0
wa_new     = 0.d0
Ufield     = 0.d0
rdiag1     = 0.d0

allocate(ds_free(ns_free_max+1))
allocate(koeff_free(ns_free_max+1))
allocate(phia_fr(numnp))
allocate(qf(numnp,2))
allocate(qf_final(numnp,ns_free_max+1))

ds_free    = 0.d0
koeff_free = 0.d0
phia_fr    = 0.d0
qf_final   = 0.d0
qf         = 0.d0

if (use_grafted.eq.1) then
    allocate(koeff_gr(ns_gr+1))
    allocate(ds_gr(ns_gr+1))
    allocate(qgr(numnp,2))

    koeff_gr   = 0.d0
    ds_gr      = 0.d0
    qgr        = 0.d0
endif

allocate(qgr_final(numnp,ns_gr+1))
qgr_final  = 0.d0

allocate(phia_gr(numnp))
phia_gr    = 0.d0

#ifdef USE_MPI
call MPI_BCAST(mumps_matrix_type, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif

!output mesh characteristics
write(iow,'(/''Mesh characteristics..'')')
write(iow,'(''   Number of mesh points (numnp):         '',I16)') numnp
write(iow,'(''   Number of elements (numel):            '',I16)') numel
write(iow,'(''   Number of nodes per element (nel):     '',I16)') nel
write(iow,'(''   Number of matrix indeces:              '',I16)') all_el

write(6  ,'(/''Mesh characteristics..'')')
write(6  ,'(''   Number of mesh points (numnp):         '',I16)') numnp
write(6  ,'(''   Number of elements (numel):            '',I16)') numel
write(6  ,'(''   Number of nodes per element (nel):     '',I16)') nel
write(6  ,'(''   Number of matrix indeces:              '',I16)') all_el

!initialize time integration scheme
call init_time(ns_free, ds_ave_free, ds_free, koeff_free)

if (use_grafted.eq.1) then
    call init_time(ns_gr, ds_ave_gr, ds_gr, koeff_gr)
endif

!initialize field
call init_field(field_in_filename, Ufield, wa)
wa_mix = wa

!initialize files in case this is not a restart
call init_files

!**************************************************************************************************************!
!                                        LOOPS FOR FIELD CONVERGENCE                                           !
!**************************************************************************************************************!
write(iow,'(/''*Initiating the simulation with '',I10,'' iterations'',/)') iterations
write(6  ,'(/''*Initiating the simulation with '',I10,'' iterations'',/)') iterations

write(iow,'(A10,1x,9(A19,1X),A16)') 'iter', 'fraction', 'adh_ten', 'n_gr_chains', 'max_error', 'std_error', &
    &                               'wa_max', 'wa_max_abs', 'wa_ave', 'wa_step'
write(6  ,'(A4, 1x,9(A14,1X),A12)') 'iter', 'fraction', 'adh_ten', 'n_gr_chains', 'max_error', 'std_error', &
    &                               'wa_max', 'wa_max_abs', 'wa_ave', 'wa_step', 'progress (%)'

t_init = get_sys_time()
iter=init_iter
do while ((iter.lt.iterations).and.(max_error.gt.max_error_tol))

    iter=iter+1

    write(iow,'(I10,1X,9(E19.9e3,1X))')              iter-1, frac, adh_ten, nch_gr, max_error, wa_std_error, &
   &                                                 wa_max, wa_max_abs, wa_ave, wa_step
    write(6  ,'(I4 ,1X,9(E14.4e3,1X))',advance='no') iter-1, frac, adh_ten, nch_gr, max_error, wa_std_error, &
   &                                                 wa_max, wa_max_abs, wa_ave, wa_step

    !flush output
    close(iow)
    open(unit=iow, file = 'scft.out.txt', position = 'append')

    !update field
    wa = wa_mix

    !perform matrix assembly
    call matrix_assemble(Rg2_per_mon_free, wa)

    !**********************SOLVE EDWARDS PDE FOR FREE CHAINS*************************!
    !initial value of propagator, qf(numnp,0) = 1.0 for all numnp
    !the initial values stored to qf_final for s=0

    do i1 = 1, numnp
       qf(i1,1)       = 1.d0
       qf_final(i1,1) = 1.d0
    enddo

    !solve
    call edwards(ds_free, ns_free_max, mumps_matrix_type, qf, qf_final)

    !*******************SOLVE EDWARDS PDE FOR GRAFTED CHAINS*************************!
    !initial value of propagator, qgr(numnp,0) = 0.0 for all numnp except for gp
    !the initial values are stored to qgr_final for s=0

    if (use_grafted.eq.1) then

        !perform matrix assembly
        call matrix_assemble(Rg2_per_mon_gr, wa)

        !re-initialize grafted propagators
        qgr = 0.d0
        qgr_final = 0.d0

        !assign initial conditions at the grafting points
        call grafted_init_cond(ns_gr, numnp, gp_filename, qgr, qgr_final)

        !solve
        call edwards(ds_gr, ns_gr, mumps_matrix_type, qgr, qgr_final)
    endif

    !*********************CONVOLUTION AND ENERGY***********************************!
    !calculate reduced segment density profiles of free and grafted chains
    call convolution(numnp, chainlen_free, ns_free, koeff_free, qf_final, qf_final, phia_fr)

    if (use_grafted.eq.1) then
        call convolution(numnp, chainlen_gr, ns_gr, koeff_gr, qgr_final, qf_final, phia_gr)
    endif

    !calculate partition function of free chains
    call part_fun(numnp, ns_free, qf_final, part_func)

    !calculated number of grafted chains
    if (use_grafted.eq.1) then
        call grafted_chains(numnp, chainlen_gr, rho_0, phia_gr, nch_gr)
    endif

    !calculate the new field
    do k1 = 1, numnp
        wa_new(k1) = kapa * (phia_fr(k1) + phia_gr(k1) - 1.d0) + Ufield(k1)
    enddo

    !calculate adhesion tension
    call adhesion_tension(qf_final, qgr_final, wa, Ufield, phia_fr, phia_gr, part_func, adh_ten)

    !*******************************************************************!
    !   COMPUTE DIFFERENCES BETWEEN OLD (wa) AND NEW (wa_new) fields    !
    !*******************************************************************!
    wa_ave       = 0.d0
    max_error    = 0.d00
    wa_std_error = 0.d00
    wa_max       = 0.d0
    wa_max_abs   = 0.d0

    do k1 = 1, numnp
        wa_ave       = wa_ave + wa_new(k1)
        max_error    = max(max_error, dabs(wa_new(k1) - wa(k1)))
        wa_std_error = wa_std_error + (wa_new(k1) - wa(k1))**2
        wa_max       = max(wa_max, wa_new(k1))
        wa_max_abs   = max(wa_max_abs, dabs(wa_new(k1)))
    enddo

    wa_std_error = SQRT(wa_std_error / float((numnp - 1)))
    wa_ave       = wa_ave / numnp

    ! TEMP NEW
    wa_ave = wa_ave * chainlen_free
    wa_max = wa_max * chainlen_free
    wa_max_abs = wa_max_abs
    max_error = max_error * chainlen_free
    wa_std_error = wa_std_error * chainlen_free

    !*******************************************************************!
    !                         APPLY MIXING RULE                         !
    !*******************************************************************!
    !original convergence scheme 1: tolis
    if (scheme_type.eq.1) then
        if (iter.eq.1) then
            frac = 1.d0
        else if (iter.eq.2) then
            frac = 1.d0 / (wa_max_abs * chainlen_free * 10.d0)
        else
            frac = min(frac * mix_coef_frac, mix_coef_kapa)
        endif

        !mixing the fields..
        do k1 = 1, numnp
            wa_mix(k1) = (1.d0 - frac) * wa(k1) + frac * wa_new(k1)
        enddo
    endif

    !experimental convergence scheme 2: 1/wa_max
    if (scheme_type.eq.2) then
        mix_tol = 0.1d0  * kapa
        frac    = frac
        wa_step = 4.d0 *exp(-dble(iter)/10.d0)

        n_outside = 0
        do k1 = 1, numnp
            if (dabs(wa_new(k1)-wa(k1)) * chainlen_free > mix_tol.and.iter.le.60) then
                n_outside = n_outside + 1

                if (wa_new(k1).gt.wa(k1)) then
                    wa_mix(k1) = wa(k1) + wa_step*rand()
                else if (wa_new(k1).lt.wa(k1)) then
                    wa_mix(k1) = wa(k1) - wa_step*rand()
                endif
            else
! APS TEMP
!                if (k1.eq.gnode_id) then
!                    wa_mix(k1) = (1.d0 - 0.1d0) * wa(k1) + 0.1d0 * wa_new(k1)
!                else
                    wa_mix(k1) = (1.d0 - frac) * wa(k1) + frac * wa_new(k1)
!                endif
            endif
        enddo

!        do k1 = 1, numnp
!            if (elem_in_q0_face(k1)) then
!                wa(k1) = -kapa
!            endif
!        enddo
    endif

    !experimental convergence scheme 3: 1/wa_max_abs
    if (scheme_type.eq.3) then
        !fraction remains constant throughout the simulation
        frac = frac
        !mixing the fields..
        do k1 = 1, numnp
            wa_mix(k1) = (1.d0 - frac) * wa(k1) + frac * wa_new(k1)
        enddo
    endif

    !*******************************************************************!
    !                           PERIODIC DUMPER                         !
    !*******************************************************************!
    !output the data every so many steps
    if (mod(iter,output_every).eq.0) then
        call periodic_dumper(qf_final, qgr_final, phia_fr, phia_gr, wa, wa_new, wa_mix)
        call export_field(wa_mix, numnp, iter)
    endif

enddo!iter

!output the data at the end of the simulation
call periodic_dumper(qf_final, qgr_final, phia_fr, phia_gr, wa, wa_new, wa_mix)
call export_field(wa_mix, numnp, iter)

write(iow,'(I10,1X,9(E19.9e3,1X))')  iter, frac, adh_ten, nch_gr, max_error, wa_std_error, &
   &                                 wa_max, wa_max_abs, wa_ave, wa_step
write(6  ,'(I4 ,1X,9(E14.4e3,1X))')  iter, frac, adh_ten, nch_gr, max_error, wa_std_error, &
   &                                 wa_max, wa_max_abs, wa_ave, wa_step

if (max_error.lt.max_error_tol) then
    write(iow,'(/''Convergence of max error'',F16.9)') max_error
    write(6  ,'(/''Convergence of max error'',F16.9)') max_error
else
    write(iow,'(/''Convergence of '',I10, '' iterations'')') iterations
    write(6  ,'(/''Convergence of '',I10, '' iterations'')') iterations
endif



write(iow,'(''-----------------------------------'')')
write(6  ,'(''-----------------------------------'')')
write(iow,'(3x,A40,E16.9)')adjl('Adhesion tension (mN/m):',40),adh_ten
write(6  ,'(3x,A40,E16.9)')adjl('Adhesion tension (mN/m):',40),adh_ten
write(iow,'(3x,A40,E16.9)')adjl('Partition function Q:'   ,40),part_func
write(6  ,'(3x,A40,E16.9)')adjl('Partition function Q:'   ,40),part_func
write(iow,'(3x,A40,E16.9)')adjl('n/n_bulk:',40)               ,nch_gr * chainlen_gr / (rho_0*volume*1.d-30)
write(6  ,'(3x,A40,E16.9)')adjl('n/n_bulk:',40)               ,nch_gr * chainlen_gr / (rho_0*volume*1.d-30)

! Please do not alter the output of the following line!
write(iow,'(3x,A40,E16.9)')adjl('grafting density (A^-2):',40),nch_gr/interf_area

t_final = get_sys_time()
write(6  ,'(3x,A40,I16)')adjl('Run duration:',40), t_final - t_init

#ifdef USE_MPI
!root will send a stop signal to the slaves
if (root) then
    flag_continue = .false.
    call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
end if

1000 call MPI_FINALIZE(ierr)
#endif
write(iow,'(/''Done!'')')
!------------------------------------------------------------------------------------------------------------------!
end program FEM_3D
