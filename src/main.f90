Program FEM_3D
!-------------------------------------------------------------------!
use xdata
use constants
use error_handing
use write_helper
use kcw
#ifdef USE_MPI
use mpistuff
#endif
!-------------------------------------------------------------------!
implicit none
#ifdef USE_MPI
include 'mpif.h'
#endif
!-------------------------------------------------------------------!
integer :: i1, k1, iter, n_outside

character(12) :: field_in_filename = 'field_in.bin'
character(40) :: field_filename_aux = ''

real(8)                              :: surf_pot
real(8)                              :: wa_max, wa_max_abs
real(8)                              :: mix_tol, wa_step, wa_ave
real(8)                              :: part_func, nch_per_area, coef
real(8)                              :: adh_ten
real(8), allocatable, dimension(:)   :: wa, wa_new, wa_mix, Ufield
real(8), allocatable, dimension(:)   :: phia_fr, phia_gr
real(8), allocatable, dimension(:,:) :: qf, qgr
real(8), allocatable, dimension(:,:) :: qf_final, qgr_final


! APS TEMP
integer :: gnode_id = 20039

!*******************************************************************!
!                           MPI SECTION                             !
!*******************************************************************!

#ifdef USE_MPI
call mpi_init(ierr)
call MPI_Comm_size ( MPI_COMM_WORLD, n_proc, ierr )
call MPI_Comm_rank ( MPI_COMM_WORLD, my_id, ierr )

if ( my_id == 0 ) then
    root = .true.
else
    root = .false.
endif

if (root) then
    write ( *, *)
    write ( *, '(a,i4,a)' ) ' MPI run with ', n_proc, ' procs.'
    write ( *, *)
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
            call mumps_sub(0,mumps_matrix_type)
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
write(iow,'(''Polumer Bulk with solid surfaces   ''/  &
          & ''---------------------------------  ''/    &
          & ''SCF theory, Gaussian string model  ''/    &
          & ''Finite Element solution with       ''/    &
          & ''successive substitutions           '')')

!Initialize the error log
ioe = 11
open(unit=ioe, file = 'error.out.txt', status='replace')
close(ioe)
!*******************************************************************!
!                       INITIALIZATION SECTION                      !
!*******************************************************************!

!parse the inputs from the datafile
call scfinout

!read the input from the mesh file and generate it
call mesh_io_3d

!allocate and initialize essential scf arrays
allocate(wa(numnp),wa_mix(numnp),wa_new(numnp),Ufield(numnp))
allocate(phia_fr(numnp),phia_gr(numnp))
allocate(rdiag1(numnp))
allocate(qf(numnp,2),qgr(numnp,2))
allocate(qf_final(numnp,ns+1),qgr_final(numnp,ns+1))     !this needs to be generalized: ns_fr vs ns_gr        

wa        = 0.d0
wa_mix    = 0.d0
wa_new    = 0.d0
Ufield    = 0.d0
phia_fr   = 0.d0
phia_gr   = 0.d0
qf        = 0.d0
qgr       = 0.d0
qf_final  = 0.d0
qgr_final = 0.d0
rdiag1    = 0.d0

!
wa_max = 0.d0
wa_ave = 0.d0
wa_step = 0.d0

#ifdef USE_MPI
call MPI_BCAST(mumps_matrix_type, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif

!*******************************************************************!
!               INITIALIZE THE TIME INTEGRATION SCHEME              !
!*******************************************************************!

!discretize the time domain and compute the simpson coeffs
allocate(ds(ns+1))
allocate(koeff(ns+1))
allocate(xs(ns+1))
ds = 0.d0
koeff = 0.d0
xs = 0.d0

!Constant step size
if (time_integration_scheme.eq.1) then
    ds = ds_ave
    call simpsonkoef_s(koeff, ds_ave, ns)
endif

!Non constant step size with quadradic interpolation
if (time_integration_scheme.eq.2) then

! APS TODO: move these inside a routine
! A)
!    ds = ds_ave
!    do k1 = 2, ns+1
!        xs(k1) = xs(k1-1) + ds_ave
!        ds(k1) = xs(k1) - xs(k1-1)
!    enddo
!
! B)

    ds(1)=0.d0
    do k1 = 2, ns+1
        xs(k1) = 0.5d0 * (1.d0 - DCOS(pi * (dble(k1)-1.d0) / dble(ns)))
        ds(k1) = xs(k1) - xs(k1-1)
    enddo

    call quadinterp_koef(koeff, xs, ds, ns)
endif

!output the mesh characteristics

write(iow,'(/''Mesh characteristics..'')')
write(iow,'(''   Number of mesh points (numnp):         '',I16)')numnp
write(iow,'(''   Number of elements (numel):            '',I16)')numel
write(iow,'(''   Number of nodes per element (nel):     '',I16)')nel
write(iow,'(''   Number of matrix indeces:              '',I16)')all_el

write(6  ,'(/''Mesh characteristics..'')')
write(6  ,'(''   Number of mesh points (numnp):         '',I16)')numnp
write(6  ,'(''   Number of elements (numel):            '',I16)')numel
write(6  ,'(''   Number of nodes per element (nel):     '',I16)')nel
write(6  ,'(''   Number of matrix indeces:              '',I16)')all_el

!*******************************************************************!
!                       INITIALIZE FIELDS                           !
!*******************************************************************!
open(unit=211, file = 'Usolid.out.txt')
do k1 = 1, numnp
    !TODO: fix distance for 3D
    distance   = xc(1,k1)
    Ufield(k1) = surf_pot(distance)
    write(211,('(E16.9,2X,E19.9)')) distance, Ufield(k1)
    if (Ufield(k1).ne.Ufield(k1)) then
        Ufield(k1) = 0.d0
        write(ERROR_MESSAGE,'(''Hamaker assumed a NaN value for x = '',E16.9,''. NaN was changed to '',E16.9)') distance, Ufield(k1)
        call exit_with_error(0,2,0,ERROR_MESSAGE)
    endif
enddo
close(211)

!*******************************************************************!
!                             READ FIELD                            !
!*******************************************************************!
if (readfield.eq.1) then
    write(iow,'(/A40,5x,A12)')adjl('*Reading field from file:',40),field_in_filename
    write(6  ,'(/A40,5x,A12)')adjl('*Reading field from file:',40),field_in_filename

    INQUIRE(FILE=field_in_filename, EXIST=FILE_EXISTS)

    if (FILE_EXISTS) then
        open(unit=655, file = field_in_filename, Form='unformatted')
    else
        write(ERROR_MESSAGE,'(''File '',A15,'' does not exist!'')')field_in_filename
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    read(655) wa
    close(655)

    !multiply field with chain length since it is divided with chain length right
    !before it is printed
    do k1 = 1, numnp
        wa(k1) = wa(k1) * chainlen
    enddo
else
    if (scheme_type.eq.2) then
        do k1 = 1, numnp
            if (elem_in_q0_face(k1)) then
                wa(k1) = -kapa
            else
                wa(k1) = 0.d0
            endif
        enddo
    else
        wa = 0.d0
    endif
endif

!Initialize the files in case this is not a restart
if (init_iter.eq.0) then
    open(unit=121, file = 'wa.out.txt', status='replace')
    close(121)
    open(unit=121, file = 'wa_new.out.txt', status='replace')
    close(121)
    open(unit=121, file = 'wa_mix.out.txt', status='replace')
    close(121)
    open(unit=121, file = 'rho.out.txt', status='replace')
    close(121)
endif

!*******************************************************************!
!                       LOOPS FOR SOLUTION                          !
!*******************************************************************!
write(iow,'(/''*Initiating the simulation with '',I10,'' iterations'',/)') iterations
write(6  ,'(/''*Initiating the simulation with '',I10,'' iterations'',/)') iterations

write(iow,'(A10,1x,8(A19,1X),A16)') 'iter', 'adh_ten', 'max_error', 'std_error', 'wa_max', 'wa_max_abs', 'wa_ave', 'wa_step', &
    &                               'fraction'
write(6  ,'(A4, 1x,8(A14,1X),A12)') 'iter', 'adh_ten', 'max_error', 'std_error', 'wa_max', 'wa_max_abs', 'wa_ave', 'wa_step', &
    &                               'fraction', 'progress (%)'

iter=init_iter
max_error=200000.

do while ((iter.lt.iterations).and.(max_error.gt.max_error_tol))
    iter=iter+1

    write(iow,'(I10,1X,8(E19.9e3,1X))')              iter-1, adh_ten, max_error, std_error, wa_max, wa_max_abs, wa_ave, &
   &                                                 wa_step, frac
    write(6  ,'(I4 ,1X,8(E14.4e3,1X))',advance='no') iter-1, adh_ten, max_error, std_error, wa_max, wa_max_abs, wa_ave, &
   &                                                 wa_step, frac

    !Flush the output
    close(iow)
    open(unit=iow, file = 'scft.out.txt', position = 'append')

    call matrix_assemble(wa)

    !**********************SOLVE EDWARDS FOR FREE CHAINS*************************!
    !Initial value of propagator, qf(numnp,0) = 1.0 for all numnp
    !The initial values stored to qf_final for s=0

    do i1 = 1, numnp
       qf(i1,1)       = 1.d0
       qf_final(i1,1) = 1.d0
    enddo

    !Solve the diffusion equation for free chains
    call edwards_film_fem(qf, qf_final)

    !*******************SOLVE EDWARDS FOR GRAFTED CHAINS*************************!
    !Initial value of propagator, qgr(numnp,0) = 0.0 for all numnp except for gp
    !The initial values are stored to qgr_final for s=0

    if (use_grafted.eq.1) then

        do i1 = 1, numnp
           qgr(i1,1)       = 0.d0
           qgr_final(i1,1) = 0.d0
        enddo

        !Specify grafting points
        if (gnode_id > numnp) then
            write(ERROR_MESSAGE,'(''ID of grafted chain ('',I10,'') is larger from numnp ('',I10,'')'')') gnode_id, numnp
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        end if

        qgr(gnode_id,1)       = 350.0d00
        qgr_final(gnode_id,1) = 350.0d00

        !Solve the diffusion equation for grafted chains
        call edwards_film_fem(qgr, qgr_final)
    endif

    !*********************CONVOLUTION AND ENERGY***********************************!  
    !Calculate reduced segment density profiles of free and grafted chains
    call convolution(qf_final, qf_final, phia_fr)

    if (use_grafted.eq.1) then
        call convolution(qgr_final, qf_final, phia_gr)
    endif

    !Calculate partition function of free chains
    call part_fun_phi(qf_final, phia_fr, part_func, nch_per_area, coef)

    !Calculate the new field
    do k1 = 1, numnp
        wa_new(k1) = kapa * (phia_fr(k1) + phia_gr(k1) - 1.d0) + Ufield(k1)
    enddo

    !Calculate adhesion tension
    call adhesion_tension(qf_final, wa, Ufield, phia_fr, phia_gr, part_func, adh_ten)

    !*******************************************************************!
    !   COMPUTE DIFFERENCES BETWEEN OLD (wa) AND NEW (wa_new) fields    !
    !*******************************************************************!

    !APS TEMP set field min equal to -kapa
    wa_ave = 0.d0
    max_error = 0.d00
    std_error = 0.d00
    wa_max=0.d0
    wa_max_abs=0.d0
    do k1 = 1, numnp
        wa_ave = wa_ave + wa_new(k1)
        max_error = max(max_error, dabs(wa_new(k1) - wa(k1)))
        std_error = std_error + (wa_new(k1) - wa(k1))**2
        wa_max = max(wa_max, wa_new(k1))
        wa_max_abs = max(wa_max_abs, dabs(wa_new(k1)))
    enddo
    std_error = SQRT(std_error / float((numnp - 1)))
    wa_ave = wa_ave / numnp
    !*******************************************************************!
    !                         APPLY MIXING RULE                         !
    !*******************************************************************!

    !   original convergence scheme 1: tolis
    if (scheme_type.eq.1) then
        if (iter.eq.1) then
            frac = 1.d0
        else if (iter.eq.2) then
            frac = 1.d0 / (wa_max_abs * 10.d0)
        else
            frac = min(frac * mix_coef_frac, mix_coef_kapa)
        endif

        !Mixing the fields..
        do k1 = 1, numnp
            wa_mix(k1) = (1.d0 - frac) * wa(k1) + frac * wa_new(k1)
        enddo
    endif

    !experimental convergence scheme 2: 1/wa_max
    if (scheme_type.eq.2) then
        mix_tol = 0.1d0  * kapa
        frac    = 0.005d0
        wa_step = 4.d0 *exp(-dble(iter)/10.d0)

        n_outside = 0
        do k1 = 1, numnp
            if (dabs(wa_new(k1)-wa(k1)) > mix_tol.and.iter.le.60) then
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
        !Fraction remains constant throughout the simulation
        frac = frac
        !Mixing the fields..
        do k1 = 1, numnp
            wa_mix(k1) = (1.d0 - frac) * wa(k1) + frac * wa_new(k1)
        enddo
    endif

    !*******************************************************************!
    !                        PERIODIC PROFILER                          !
    !*******************************************************************!

    !Output the data every so many steps
    if (mod(iter,output_every).eq.0) then
        open (unit=121, file = 'wa.out.txt', status='unknown', position='append')
        open (unit=122, file = 'wa_new.out.txt', status='unknown', position='append')
        open (unit=123, file = 'wa_mix.out.txt', status='unknown', position='append')
        open (unit=124, file = 'rho_free.out.txt', status='unknown', position='append')
        open (unit=125, file = 'rho_grafted.out.txt', status='unknown', position='append')
        do k1 = 1, numnp, print_ev
            write(121,'(E19.9e3)',advance='no') wa(k1)
            write(122,'(E19.9e3)',advance='no') wa_new(k1)
            write(123,'(E19.9e3)',advance='no') wa_mix(k1)
            write(124,'(E19.9e3)',advance='no') phia_fr(k1)
            write(125,'(E19.9e3)',advance='no') phia_gr(k1)
        enddo
        write(121,*)
        write(122,*)
        write(123,*)
        write(124,*)
        write(125,*)
        close(121)
        close(122)
        close(123)
        close(124)
        close(125)

        open (unit=120, file = 'rho_reduced.out.txt')
        write(120,'(A13,8(A19))') 'np','x','y','z','phi_f','phi_g','wa','wa_new','wa_mix'
        do k1 = 1, numnp
            write(120,'(I13,8(E19.9e3))') k1, xc(1,k1), xc(2,k1), xc(3,k1), phia_fr(k1), &
   &                                      phia_gr(k1), wa(k1), wa_new(k1), wa_mix(k1)
        enddo
        close(120)

        ! Print the restricted partition functions
        call qprint(qf_final,"free")

        if (use_grafted.eq.1) then
            call qprint(qgr_final,"graf")
        endif

    endif

    wa = wa_mix

    !*******************************************************************!
    !               EXPORT FIELD TO BINARY FILE FOR RESTART             !
    !*******************************************************************!

    if (mod(iter,output_every).eq.0) then
        !Normalize the field by dividing it with chain length
        do k1 = 1, numnp
            wa_mix(k1) = wa_mix(k1) / chainlen
        enddo

        field_filename_aux = ""
        write(field_filename_aux,'(''field_'',I4.4,''.out.bin'')')iter
        open(unit=655, file = field_filename_aux, Form='unformatted')
        write(655) wa_mix
        close(655)

        field_filename_aux = ""
        write(field_filename_aux,'(''field_comsol_'',I4.4,''.out.txt'')')iter
        open(unit=123, file = field_filename_aux)

        do k1 = 1, numnp
            write(123,'(E16.9,2(2X,E16.9),(2X,E19.9e3))') xc(1,k1), xc(2,k1), xc(3,k1), -wa_mix(k1)
        enddo

        close(123)
    endif
enddo!iter


write(iow,'(I10,1X,8(E19.9e3,1X))')              iter, adh_ten, max_error, std_error, wa_max, wa_max_abs, wa_ave, &
   &                                                 wa_step, frac
write(6  ,'(I4 ,1X,8(E14.4e3,1X))',advance='no') iter, adh_ten, max_error, std_error, wa_max, wa_max_abs, wa_ave, &
   &                                                 wa_step, frac

!write(iow,'(I10,1X,6(E19.9e3,1X))')              iter, adh_ten, max_error, std_error, wa_max, wa_max_abs, frac
!write(6  ,'(I4 ,1X,6(E14.4e3,1X))',advance='no') iter, adh_ten, max_error, std_error, wa_max, wa_max_abs, frac

if (max_error.lt.max_error_tol) then
    write(iow,'(/''Convergence of max error'',F16.9)') max_error
else
    write(iow,'(/''Convergence of '',I10, '' iterations'')') iterations
endif

write(iow,'(''-----------------------------------'')')
write(iow,'(''Adhesion tension (mN/m) '',E16.9)') adh_ten
write(iow,'(''Partition function Q =  '',E16.9)') part_func
write(iow,'(''            n/n_bulk =  '',E16.9)') nch_per_area * chainlen / (rho_0*volume*1.d-30)

#ifdef USE_MPI
!Root will send a stop signal to the slaves
if (root) then
    flag_continue = .false.
    call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
end if

1000 call MPI_FINALIZE(ierr)
#endif
write(*,'(/''Done!'')')
!------------------------------------------------------------------------------------------------------------------!
end program FEM_3D
