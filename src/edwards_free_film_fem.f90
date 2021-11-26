subroutine edwards_free_film_fem
!----------------------------------------------------------------------------------------------------------!
use xdata
use constants
use mdata
use kcw
#ifdef USE_MPI
use mpistuff
#endif
!----------------------------------------------------------------------------------------------------------!
implicit none
#ifdef USE_MPI
include 'mpif.h'
#endif
!----------------------------------------------------------------------------------------------------------!
integer :: i, j, f, idummy, time_step

logical, dimension(numnp) :: r_true
logical, dimension(nel*numel) :: set_diag_to_one

#ifdef PRINT_AFULL
real(8), allocatable, dimension(:,:) :: A_full
#endif

!************************INITIAL CONDITIONS*************************!
r_true = .false.

!Initial value of propagator, q(x,0) = 1.0 for all x
!The initial values stored to qf_final for s=0

do i1 = 1,numnp
    qf(i1,1) = 1.d0
    qf_final(i1,1) = 1.d0
enddo

g_m%value  = c_m%value + ds*(k_m%value + w_m%value)

rh_m%value = c_m%value

!************************BOUNDARY CONDITIONS************************!
do j = 1, fcel
    do i1 = 1, n_dirichlet_faces
        if (fcentity(j)==ids_dirichlet_faces(i1))then
            do i = 1, fcnum
                idummy= fcelement(i,j)
                r_true(idummy) = .True.
            enddo
        endif
    enddo
enddo

set_diag_to_one=.true.

do i1 = 1, all_el
    f = g_m%col(i1)
    i = g_m%row(i1)
#if defined(MSYMDEFPOS) || defined(MSYMGEN)
    if (i > f) then
        g_m%value(i1) = 0.d0
    endif
    if (r_true(i).or.r_true(f)) then
#else
    if (r_true(i)) then
#endif
        g_m%value(i1) = 0.d0
        if (i==f.and.set_diag_to_one(i)) then
            g_m%value(i1) = 1.d0
            set_diag_to_one(i)=.false.
        endif
    endif
enddo

!***********************DETERMINE NON-ZERO ENTRIES*******************!
NNZ = 0
do i = 1, all_el
    if (abs(g_m%value(i)) > tol) then
        NNZ = NNZ + 1
    endif
enddo

allocate(A_m%value(NNZ))
allocate(A_m%col(NNZ))
allocate(A_m%row(NNZ))

NNZ = 0
do i = 1, all_el
    if (abs(g_m%value(i)) > tol) then
        NNZ = NNZ + 1

        A_m%value(NNZ) = g_m%value(i)
        A_m%row(NNZ)   = g_m%row(i)
        A_m%col(NNZ)   = g_m%col(i)
    endif
enddo

#ifdef PRINT_AFULL
allocate(A_full(numnp,numnp))
A_full = 0.d0
do i1 = 1, NNZ
   f = A_m%col(i1)
   i = A_m%row(i1)
   A_full(i,f)=A_m%value(i1)
enddo
open(unit=255, file = 'A_full.out.txt')
do i = 1, numnp
   write(255,*)(A_full(i,f), f = 1, numnp)
enddo
close(255)
deallocate(A_full)
#endif


!************************START TRANSIENT SOLUTION*********************! 
write(6,'(A11)' ,advance='no') 'time_step = '
do time_step = 2, ns+1

#ifdef USE_MPI
    ! Send a continue (.true.) signal to the slaves
    call MPI_BCAST(.true., 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    write(6,'(I3,A2)', advance='no') time_step, ", "
    !*******************************************************************!
    !                   FINITE ELEMENT METHOD                           !
    !*******************************************************************!
    rdiag1 = 0.

    do i1 = 1, all_el
        i = rh_m%row(i1)
        j = rh_m%col(i1)

        rdiag1(i) = rdiag1(i) + rh_m%value(i1)*qf(j,1)
    enddo

    do i = 1, numnp
        if (r_true(i)) rdiag1(i) = 0.
    enddo

    call mumps_sub(numnp)

    do i1 = 1,numnp
         qf(i1,2) = rdiag1(i1)
    enddo

    !save propagators for convolution
    do i1 = 1,numnp
        qf_final(i1,time_step) = qf(i1,2)
        qf(i1,1) = qf(i1,2)
    enddo

enddo !time_step
write(6,*)

return
!----------------------------------------------------------------------------------------------------------!	
end subroutine edwards_free_film_fem
