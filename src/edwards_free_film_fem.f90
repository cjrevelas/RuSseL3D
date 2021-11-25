subroutine edwards_free_film_fem
!----------------------------------------------------------------------------------------------------------!
use xdata
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
integer :: i, j, f, idummy!, fc

real(8) :: t_0, t_1, t_2, t_3, t_4

logical, dimension(numnp) :: r_true
!times
write(54, '(A9,A9,A9,A18  ,A18  ,A18  ,A18  )')'time','NZ','N','BC+G+R',   'NZcalc',   'MUMPS',    'cp_sol'
write(54, '(A9,A9,A9,A9,A9,A9,A9,A9,A9,A9,A9)')''    , '' , '','min','sec','min','sec','min','sec','min','sec'

!************************INITIAL CONDITIONS*************************!
call CPU_TIME(t_0)

r_true = .false.

!Initial value of propagator, q(x,0) = 1.0 for all x
!The initial values stored to qf_final for s=0

do i1 = 1,numnp
    qf(i1,1) = 1.d0
    qf_final(i1,1) = 1.d0
enddo

g_m%value  = c_m%value + ds*(k_m%value + w_m%value)

rh_m%value = c_m%value

call CPU_TIME(t_1)

!************************BOUNDARY CONDITIONS************************!
do j = 1, fcel
    if ((fcentity(j)==3).or.(fcentity(j)==4))then
        do i = 1, fcnum
            idummy= fcelement(i,j)
            r_true(idummy) = .True.
        enddo
    endif
enddo

! APS 16/08/19: OPTIMIZE

! new section
do i1 = 1, all_el
    f = g_m%col(i1)
    i = g_m%row(i1)
    if (r_true(i)) then
        g_m%value(i1) = 0.
        if (i==f) then
            g_m%value(i1) = 1.
        endif
    endif
enddo
!/ new section

! old section
!do i = 1, numnp
!    if (r_true(i)) then
!        do f = 1, numnp
!           do i1 = 1, all_el
!               if (f==rh_m%col(i1).and.i==rh_m%row(i1)) then
!                    g_m%value(i1) = 0.
!                    if (i==rh_m%col(i1).and. i==rh_m%row(i1)) then
!                        g_m%value(i1) = 1.
!                    endif
!               endif
!            enddo
!        enddo
!    endif
!enddo
!/ old section

call CPU_TIME(t_2)
!***********************DETERMINE NON-ZERO ENTRIES*******************!
NNZ = 0
do i = 1, all_el
    !APS 16/08/19: I set abs(g_m%value(i))>1.e-8 instead of g_m%value(i)/=0.
    if (abs(g_m%value(i))>1.e-8) then
        NNZ = NNZ + 1
    endif
enddo

allocate(A_m%value(NNZ))
allocate(A_m%col(NNZ))
allocate(A_m%row(NNZ))

NNZ = 0
do i = 1, all_el
    if (abs(g_m%value(i))>1.e-8) then
        NNZ = NNZ + 1

        A_m%value(NNZ) = g_m%value(i)
        A_m%row(NNZ)   = g_m%row(i)
        A_m%col(NNZ)   = g_m%col(i)
    endif
enddo

call CPU_TIME(t_3)

!************************START TRANSIENT SOLUTION*********************! 
do time_step = 2, ns+1

#ifdef USE_MPI
    ! Send a continue (.true.) signal to the slaves
    call MPI_BCAST(.true., 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    write(6,*) 'time_step =', time_step      
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

#ifdef USE_MPI
! Send a stop (.false.) signal to the slaves
call MPI_BCAST(.false., 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

call CPU_TIME(t_4)

!*********************************************************************! 
write(54, '(i9,                 i9,              i9,                                 &
    &       i9,                 f9.1,            i9,              f9.1,              &
    &       i9,                 f9.1,            i9,              f9.1)')            &
    &       time_step,          NNZ,             numnp,                              &
    &       int(t_1-t_0)/60, mod((t_1-t_0),60.), int(t_2-t_1)/60, mod((t_2-t_1),60.),&
    &       int(t_3-t_2)/60, mod((t_3-t_2),60.), int(t_4-t_3)/60, mod((t_4-t_3),60.)

return
!----------------------------------------------------------------------------------------------------------!	
end subroutine edwards_free_film_fem
