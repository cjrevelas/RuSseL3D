Program FEM_3D
!-------------------------------------------------------------------!
use xdata
use constants
use mdata
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
real(8) :: surf_pot
real(8) :: wa_max
logical :: zero_field
integer :: print_ev
!-------------------------------------------------------------------!

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
!
! the slaves will enter the mumps subroutine until they receive a stop
! signal from master proc
if (.not.root) then
    do while (.true.)
        call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        if (flag_continue) then
            call mumps_sub(0)
        else
            exit
        endif
    end do
    goto 1000
endif
#endif
!*******************************************************************!
!                       INITIALIZATION SECTION                      !
!*******************************************************************!
ior = 55
iow = 10

!TODO ADD THESE TO INITIAL LOG
write(6,*)"Dirichel BCs have been applied in faces # and #"

#if defined(MSYMGEN)
write(6,*)"The matrices are treated as GENERAL SYMMETRIC"
#endif
#if defined(MSYMDEFPOS)
write(6,*)"The matrices are treated as SYMMETRIC POSITIVE DEFINITE"
#endif
#if !defined(MSYMGEN) && !defined(MSYMDEFPOS)
write(6,*)"The matrices are treated as UNSYMMETRIC"
#endif

write(6,*) "Initializing mesh.."
call mesh_io_3d

write(6,*) "Reading gaussdat parameter values.."
call scfinout

! Initialize the Simpson coefficients
call simpsonkoef_s

!*******************************************************************!
!                       INITIALIZE FIELDS                           !
!*******************************************************************!
wa = 0.d0

open(unit=211, file = 'Usolid.out.txt')
do k1 = 1, numnp
    !TODO: fix distance for 3D
    distance   = xc(1,k1)
    Ufield(k1) = surf_pot(distance)
    write(211,('(E16.9,2X,E19.9)')) distance, Ufield(k1)
enddo
close(211)

zero_field = .true.
if (show.eq.1) then
    write(*,*)"Field will be read from a file!"
    open(unit=21, file = 'field.in.bin', Form='unformatted')
    read(21) wa
    close(21)
#ifdef REDUCE_W_CHLEN
    do k1 = 1, numnp
        wa(k1) = wa(k1) * chainlen
    enddo
#endif
    zero_field = .false.
endif
!*******************************************************************!
!                       LOOPS FOR SOLUTION                          !
!*******************************************************************!
write(6,*) 'Iteration,  Adh. tension (mN/m),  error(beta N w)'
write(iow,'(A10,2X,A16,2X,A16)') 'iter', 'adh_ten', 'error'

iter=0
error=200000.

do while ((iter.lt.iterations).and.(error.gt.max_error))
    iter=iter+1

    call matrix_assemble

    ! Solve the diffusion equation
    call edwards_free_film_fem

    ! Calculate the partition function and reduced density
    call part_fun_phi

    ! Calculate the new field
    do k1 = 1, numnp
        wa_new(k1) = kapa * (phia_new(k1)- 1.d0) + Ufield(k1)
    enddo

    ! Calculate adhesion tension
    call adhesion_tension

    !*******************************************************************!
    !                     PERIODIC EXPORT OF PROFILES                   !
    !*******************************************************************!

!APS TEMP Generate the headers of files
    print_ev = 3
    if (zero_field.and.iter.eq.1) then
        open (unit=120, file = 'rho.out.txt')
        do k1 = 1, numnp, print_ev
            write(120,'(E13.5)',advance='no') xc(3,k1)
        enddo
        close(120)
        open (unit=120, file = 'wa.out.txt')
        do k1 = 1, numnp, print_ev
            write(120,'(E13.5)',advance='no') xc(3,k1)
        enddo
        close(120)
        open (unit=120, file = 'wa_new.out.txt')
        do k1 = 1, numnp, print_ev
            write(120,'(E13.5)',advance='no') xc(3,k1)
        enddo
        close(120)
        open (unit=120, file = 'wa_mix.out.txt')
        do k1 = 1, numnp, print_ev
            write(120,'(E13.5)',advance='no') xc(3,k1)
        enddo
        close(120)
    endif
!/APS TEMP SECTION

    if (pr_on==1) then
        call qprint
    endif

    !*******************************************************************!
    !                     APPLY CONVERGENCE CRITERION                   !
    !*******************************************************************!
    error = 0.d00
    do k1 = 1, numnp
        error = error + (dabs(wa_new(k1) - wa(k1)))
    enddo
    error = error / numnp
    !*******************************************************************!
    !                         APPLY MIXING RULE                         !
    !*******************************************************************!

    wa_max=0.d0
    ! APS TEMP dabs(wa_new) might not work properly for wa_max < kapa!!!
    do k1 = 1, numnp
        wa_max = max(wa_max, dabs(wa_new(k1)))
    enddo

    !   original convergence scheme
    if (scheme_type.eq.1) then
        if (.not.zero_field) then
            fraction = min(fraction * mix_coef_frac, 1.d0 / kapa / mix_coef_kapa)
        endif
        if (zero_field.and.iter.eq.1) then
            fraction = 1.d0
        endif
        if (zero_field.and.iter.eq.2) then
            fraction = 1.d0 / (wa_max * 100.d0)
            zero_field = .false.
        endif
    endif

    ! experimental scheme
    if (scheme_type.eq.2) then
        fraction = min(1.d0 / wa_max / mix_coef_frac, 1.d0 / kapa / mix_coef_kapa)

        if (zero_field) then
            fraction = 1.d0
            zero_field = .false.
        endif
    endif

    do k1 = 1, numnp
        wa_mix(k1) = (1.d0-fraction) * wa(k1) + fraction * wa_new(k1)
    enddo

    write(6  ,'(I10,2X,E16.9,2X,E16.9,2X,E16.9,2X,E16.9)') iter, adh_ten, error, wa_max, fraction
    write(iow,'(I10,2X,E16.9,2X,E16.9,2X,E16.9,2X,E16.9)') iter, adh_ten, error, wa_max, fraction


    if (mod(iter,1).eq.0) then
        open (unit=121, file = 'wa.out.txt', status='unknown', position='append')
        open (unit=122, file = 'wa_new.out.txt', status='unknown', position='append')
        open (unit=123, file = 'wa_mix.out.txt', status='unknown', position='append')
        open (unit=124, file = 'rho.out.txt', status='unknown', position='append')
        do k1 = 1, numnp, print_ev
            write(121,'(E13.5)',advance='no') wa(k1)
            write(122,'(E13.5)',advance='no') wa_new(k1)
            write(123,'(E13.5)',advance='no') wa_mix(k1)
            write(124,'(E13.5)',advance='no') phia_new(k1)
        enddo
        write(121,*)
        write(122,*)
        write(123,*)
        write(124,*)
        close(121)
        close(122)
        close(123)
        close(124)

        open (unit=120, file = 'rho_reduced.out.txt')
        write(120,'(8(A13))') 'np','x','y','z','phi','wa','wa_new','wa_mix'
        do k1 = 1, numnp
            write(120,'(I13,7(E13.5))') k1, xc(1,k1), xc(2,k1), xc(3,k1), phia_new(k1), &
   &                                    wa(k1), wa_new(k1), wa_mix(k1)
        enddo
        close(120)
    endif

    wa = wa_mix

    ! Print the field for restart
#ifdef REDUCE_W_CHLEN
    do k1 = 1, numnp
        wa_mix(k1) = wa_mix(k1) / chainlen
    enddo
#endif
    open(unit=21, file = 'field.out.bin', Form='unformatted')
    write(21) wa_mix
    close(21)

enddo!iter

if (error.lt.max_error) then
    write(iow,'(/''Convergence of max error'',F16.9)') error
else
    write(iow,'(/''Convergence of '',I10, '' iterations'')') iterations
endif

write(iow,'(''-----------------------------------'')')
write(iow,'(''Adhesion tension (mN/m) '',E16.9)') adh_ten
write(iow,'(''Partition function Q =  '',E16.9)') part_func
write(iow,'(''            n/n_bulk =  '',E16.9)') nch_per_area * chainlen / (rho_0*volume*1.d-30)

#ifdef USE_MPI
! Root will send a stop signal to the slaves
if (root) then
    flag_continue = .false.
    call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
end if

1000 call MPI_FINALIZE(ierr)
#endif
write(*,*)"Done!"
!------------------------------------------------------------------------------------------------------------------!
end program FEM_3D
