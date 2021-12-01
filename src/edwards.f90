subroutine edwards(ds, ns, mumps_matrix_type, q, q_final)
!----------------------------------------------------------------------------------------------------------!
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
integer, intent(in) :: ns, mumps_matrix_type
integer             :: i, j, i1, time_step
integer             :: get_sys_time, t_init, t_final

real(8), intent(in), dimension(ns+1)          :: ds
real(8), intent(inout), dimension(numnp,2)    :: q
real(8), intent(inout), dimension(numnp,ns+1) :: q_final
!----------------------------------------------------------------------------------------------------------!
call dirichlet(ds, ns, mumps_matrix_type)

t_init = get_sys_time()

!************************START TRANSIENT SOLUTION*********************!
do time_step = 2, ns+1

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

    call mumps_sub(mumps_matrix_type)

    do i1 = 1,numnp
         q(i1,2) = rdiag1(i1)
    enddo

    !save propagators for convolution
    do i1 = 1,numnp
        q_final(i1,time_step) = q(i1,2)
        q(i1,1) = q(i1,2)
    enddo
enddo !time_step


deallocate(A_m%value)
deallocate(A_m%col)
deallocate(A_m%row)

write(6,'(1x,I3)',advance='no') 100

t_final = get_sys_time()

write(6,'('' :'',I6)') t_final - t_init

return
!----------------------------------------------------------------------------------------------------------!
end subroutine edwards
