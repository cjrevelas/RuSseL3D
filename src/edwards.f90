subroutine edwards(ds, q, q_final)
!----------------------------------------------------------------------------------------------------------!
use xdata
use constants
use kcw
#ifdef USE_MPI
use mpistuff
#endif
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
#ifdef USE_MPI
include 'mpif.h'
#endif
!----------------------------------------------------------------------------------------------------------!
integer :: i, j, f, i1, time_step

logical, dimension(nel*numel) :: set_diag_to_one

real(8), intent(in), dimension(ns+1)          :: ds
real(8), intent(inout), dimension(numnp,2)    :: q
real(8), intent(inout), dimension(numnp,ns+1) :: q_final

#ifdef PRINT_AFULL
real(8), allocatable, dimension(:,:) :: A_full
#endif
!----------------------------------------------------------------------------------------------------------!

#ifdef VARIABLE_DS_SCHEME
do time_step = 2, ns+1
    F_m%g  = F_m%c + ds(time_step)*(F_m%k + F_m%w)
#else
    F_m%g  = F_m%c + ds(1)*(F_m%k + F_m%w)
#endif

    F_m%rh = F_m%c


!APS TODO: The following sections must be removed from this routine
!          and placed within matrix_assembly.f90

!
! ---START OF SECTION TO BE MOVED---
!
!

    set_diag_to_one=.true.

    !in case the matrix is symmetric remove the zero lines and rows
    !diagonal componets with Dirichlet BC q=0.
    if (mumps_matrix_type.eq.1.or.mumps_matrix_type.eq.2) then
        do i1 = 1, all_el
            f = F_m%col(i1)
            i = F_m%row(i1)
            if (i > f) then
                F_m%g(i1) = 0.d0
            endif
            if (elem_in_q0_face(i).or.elem_in_q0_face(f)) then
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
            if (elem_in_q0_face(i)) then
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
!
!
! ---END OF SECTION TO BE MOVED---
!
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
    !APS TEMP
#ifndef VARIABLE_DS_SCHEME
    do time_step = 2, ns+1
#endif

#ifdef USE_MPI
        !send a continue (.true.) signal to the slaves
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

            rdiag1(i) = rdiag1(i) + F_m%rh(i1)*q(j,1)
        enddo

        do i = 1, numnp
            if (elem_in_q0_face(i)) rdiag1(i) = 0.
        enddo

        call mumps_sub(numnp, mumps_matrix_type)

        do i1 = 1,numnp
             q(i1,2) = rdiag1(i1)
        enddo

        !save propagators for convolution
        do i1 = 1,numnp
            q_final(i1,time_step) = q(i1,2)
            q(i1,1) = q(i1,2)
        enddo



    !APS TEMP - The A_m matrices must be deallocated to prevent memory leak

#ifdef VARIABLE_DS_SCHEME
    deallocate(A_m%value)
    deallocate(A_m%col)
    deallocate(A_m%row)
#endif

enddo !time_step

write(6,'(1x,I3)') 100

return
!----------------------------------------------------------------------------------------------------------!
end subroutine edwards
