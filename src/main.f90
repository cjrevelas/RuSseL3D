Program FEM_3D
!-------------------------------------------------------------------!
use xdata
use constants
use error_handing
use write_helper
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
character(12) :: field_filename = 'field.in.bin'
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
    ! receive the matrix type from root
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

! Initialize the error log
ioe = 11
open(unit=ioe, file = 'error.out.txt', status='replace')
close(ioe)
!*******************************************************************!
!                       INITIALIZATION SECTION                      !
!*******************************************************************!

call scfinout

#ifdef USE_MPI
call MPI_BCAST(mumps_matrix_type, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif

call mesh_io_3d

call simpsonkoef_s

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

zero_field = .true.
if (readfield.eq.1) then
    write(iow,'(/A40,5x,A12)')adjl('*Reading field from file:',40),field_filename
    write(6  ,'(/A40,5x,A12)')adjl('*Reading field from file:',40),field_filename

    INQUIRE(FILE=field_filename, EXIST=FILE_EXISTS)

    if (FILE_EXISTS) then
        open(unit=655, file = field_filename, Form='unformatted')
    else
        write(ERROR_MESSAGE,'(''File '',A15,'' does not exist!'')')field_filename
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    read(655) wa
    close(655)
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
write(iow,'(/''*Initiating the simulation with '',I10,'' iterations'',/)') iterations
write(6  ,'(/''*Initiating the simulation with '',I10,'' iterations'',/)') iterations

write(iow,'(A10,1x,5(A16,1X),A16)') 'iter', 'adh_ten', 'max_error', 'std_error', 'wa_max', 'fraction'
write(6  ,'(A4, 1x,5(A11,1X),A12)') 'iter', 'adh_ten', 'max_error', 'std_error', 'wa_max', 'fraction', 'progress (%)'

iter=0
max_error=200000.

do while ((iter.lt.iterations).and.(max_error.gt.max_error_tol))
    iter=iter+1

    write(iow,'(I10,1X,5(E16.9,1X))')              iter-1, adh_ten, max_error, std_error, wa_max, fraction
    write(6  ,'(I4 ,1X,5(E11.4,1X))',advance='no') iter-1, adh_ten, max_error, std_error, wa_max, fraction

    call matrix_assemble

    ! Solve the diffusion equation
    call edwards_free_film_fem

    ! Calculate the partition function and reduced density
    call part_fun_phi

    ! Calculate the new field
    do k1 = 1, numnp
        wa_new(k1) = kapa * (phia_new(k1) - 1.d0) + Ufield(k1)
    enddo

    ! Calculate adhesion tension
    call adhesion_tension

    !*******************************************************************!
    !   COMPUTE DIFFERENCES BETWEEN OLD (wa) AND NEW (wa_new) fields    !
    !*******************************************************************!
    max_error = 0.d00
    do k1 = 1, numnp
        max_error = max(max_error, dabs(wa_new(k1) - wa(k1)))
    enddo

    std_error = 0.d00
    do k1 = 1, numnp
        std_error = std_error + (wa_new(k1) - wa(k1))**2
    enddo
    std_error = SQRT(std_error / float((numnp - 1)))

    wa_max=0.d0
    do k1 = 1, numnp
        wa_max = max(wa_max, dabs(wa_new(k1)))
    enddo

    !*******************************************************************!
    !                         APPLY MIXING RULE                         !
    !*******************************************************************!

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

    !*******************************************************************!
    !                        PERIODIC PROFILER                          !
    !*******************************************************************!

    ! Initialize the files in case this is not a restart
    if (zero_field) then
        open(unit=121, file = 'wa.out.txt', status='replace')
        close(121)
        open(unit=121, file = 'wa_new.out.txt', status='replace')
        close(121)
        open(unit=121, file = 'wa_mix.out.txt', status='replace')
        close(121)
        open(unit=121, file = 'rho.out.txt', status='replace')
        close(121)
    endif

    ! Output the data every so many steps
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

        call qprint

    endif

    wa = wa_mix

    !*******************************************************************!
    !               EXPORT FIELD TO BINARY FILE FOR RESTART             !
    !*******************************************************************!
#ifdef REDUCE_W_CHLEN
    do k1 = 1, numnp
        wa_mix(k1) = wa_mix(k1) / chainlen
    enddo
#endif
    open(unit=655, file = 'field.out.bin', Form='unformatted')
    write(655) wa_mix
    close(655)

    ! Flush the output
    close(iow)
    open(unit=iow, file = 'scft.out.txt', position = 'append')

enddo!iter


write(iow,'(I10,1X,5(E16.9,1X))') iter, adh_ten, max_error, std_error, wa_max, fraction
write(6  ,'(I4 ,1X,5(E11.4,1X))',advance='no') iter, adh_ten, max_error, std_error, wa_max, fraction

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
! Root will send a stop signal to the slaves
if (root) then
    flag_continue = .false.
    call MPI_BCAST(flag_continue, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
end if

1000 call MPI_FINALIZE(ierr)
#endif
write(*,'(/''Done!'')')
!------------------------------------------------------------------------------------------------------------------!
end program FEM_3D
