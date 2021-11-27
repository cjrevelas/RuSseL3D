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

F_m%g  = F_m%c + ds*(F_m%k + F_m%w)

F_m%rh = F_m%c

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

! In case the matrix is symmetric remove the zero the lines and rows
! diagonal componets with Dirichlet BC q=0.
if (mumps_matrix_type.eq.1.or.mumps_matrix_type.eq.2) then
    do i1 = 1, all_el
        f = F_m%col(i1)
        i = F_m%row(i1)
        if (i > f) then
            F_m%g(i1) = 0.d0
        endif
        if (r_true(i).or.r_true(f)) then
            F_m%g(i1) = 0.d0
            if (i==f.and.set_diag_to_one(i)) then
                F_m%g(i1) = 1.d0
                set_diag_to_one(i)=.false.
            endif
        endif
    enddo
endif

if (mumps_matrix_type.eq.0) then
    do i1 = 1, all_el
        f = F_m%col(i1)
        i = F_m%row(i1)
        if (r_true(i)) then
            F_m%g(i1) = 0.d0
            if (i==f.and.set_diag_to_one(i)) then
                F_m%g(i1) = 1.d0
                set_diag_to_one(i)=.false.
            endif
        endif
    enddo
endif

!***********************DETERMINE NON-ZERO ENTRIES*******************!
NNZ = 0
do i = 1, all_el
    if (abs(F_m%g(i)) > tol) then
        NNZ = NNZ + 1
    endif
enddo

allocate(A_m%value(NNZ))
allocate(A_m%col(NNZ))
allocate(A_m%row(NNZ))

NNZ = 0
do i = 1, all_el
    if (abs(F_m%g(i)) > tol) then
        NNZ = NNZ + 1

        A_m%value(NNZ) = F_m%g(i)
        A_m%row(NNZ)   = F_m%row(i)
        A_m%col(NNZ)   = F_m%col(i)
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

do time_step = 2, ns+1

#ifdef USE_MPI
    ! Send a continue (.true.) signal to the slaves
    call MPI_BCAST(.true., 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    if (ns.ge.10.and.mod(time_step+1,ns/10).eq.0) then
        write(6,'(I3,1x)',advance='no') nint((time_step-2.d0)/ns*100.d0)
    elseif (ns.lt.10.and.mod(time_step+1,1).eq.0) then
        write(6,'(I3,1x)',advance='no') nint((time_step-2.d0)/ns*100.d0)
    endif
    !*******************************************************************!
    !                   FINITE ELEMENT METHOD                           !
    !*******************************************************************!
    rdiag1 = 0.

    do i1 = 1, all_el
        i = F_m%row(i1)
        j = F_m%col(i1)

        rdiag1(i) = rdiag1(i) + F_m%rh(i1)*qf(j,1)
    enddo

    do i = 1, numnp
        if (r_true(i)) rdiag1(i) = 0.
    enddo

    call mumps_sub(numnp, mumps_matrix_type)

    do i1 = 1,numnp
         qf(i1,2) = rdiag1(i1)
    enddo

    !save propagators for convolution
    do i1 = 1,numnp
        qf_final(i1,time_step) = qf(i1,2)
        qf(i1,1) = qf(i1,2)
    enddo

enddo !time_step

write(6,'(1x,I3)') 100

return
!----------------------------------------------------------------------------------------------------------!	
end subroutine edwards_free_film_fem
