subroutine mumps_sub(numnp)
!This file is part of MUMPS 5.2.1, released on Fri Jun 14 14:46:05 UTC2019
!--------------------------------------------------------------------------!
use kcw
use xdata
use mpistuff
!--------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------!
include 'mpif.h'
include 'dmumps_struc.h'
!--------------------------------------------------------------------------!
type(DMUMPS_STRUC) :: mumps_par

integer :: i, i8, i9, numnp

!--------------------------------------------------------------------------!

!Define a communicator for the package
mumps_par%COMM = MPI_COMM_WORLD

!Initialize an instance of the package for LU-factorization
mumps_par%PAR  = 1  !working host processor
mumps_par%SYM  = 0  !general case: non-symmetric left-hand-side matrix
mumps_par%JOB  = -1

call DMUMPS(mumps_par)

if (mumps_par%INFOG(1).lt.0) then
    write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",                            &
                             "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                             "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
    goto 500
endif

mumps_par%ICNTL(1) = -1
mumps_par%ICNTL(2) = -1
mumps_par%ICNTL(3) = -1
mumps_par%ICNTL(4) = -1
!mumps_par%ICNTL(13) = 1
!mumps_par%CNTL(1) = 0


!Define problem on the host processor(id = 0)
if (mumps_par%MYID.eq.0) then
    !read(5,*) mumps_par%N
    !read(5,*) mumps_par%NNZ
    mumps_par%N   = numnp
    mumps_par%NNZ = non_zero

    allocate(mumps_par%IRN( mumps_par%NNZ))
    allocate(mumps_par%JCN(mumps_par%NNZ))
    allocate(mumps_par%A(mumps_par%NNZ))
    allocate(mumps_par%RHS(mumps_par%N))
    !allocate(u1(mumps_par%N))

    mumps_par%IRN = A_m%row
    mumps_par%JCN = A_m%col
    mumps_par%A   = A_m%value
   
    do i = 1, mumps_par%N
        mumps_par%RHS(i) = rdiag1(i)
    enddo
endif
   
!DEBUG       
!if (time_step==2) then
    !open(unit=55, file = 'mumpsinside.txt')

    !do i8 = 1, mumps_par%NNZ
       !write(55,*) mumps_par%IRN(i8), mumps_par%JCN(i8), mumps_par%A(i8)
    !enddo

    !do i8 = 1, mumps_par%N
    !   write(55,*) mumps_par%RHS(i8)
    !enddo 

    !close(55)
!endif 
!DEBUG

!Call package for matrix factorization and solution
mumps_par%JOB = 6

call DMUMPS(mumps_par)
 
if (mumps_par%INFOG(1).lt.0) then
    write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",                            &
                             "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                             "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
    goto 500
endif
   
!Solution has been assembled on the host
if (mumps_par%MYID.eq.0 ) then
    !open(unit=53, file = 'u1.txt')
    do i = 1, mumps_par%N
        rdiag1(i) = mumps_par%RHS(i)
        !write(53,*) i, u1(i)
    enddo
endif

!Deallocate user data
 if (mumps_par%MYID.eq.0 ) then
    deallocate(mumps_par%IRN)
    deallocate(mumps_par%JCN)
    deallocate(mumps_par%A)
    deallocate(mumps_par%RHS)
endif

!Destruct the instance (deallocate internal data structures)
mumps_par%JOB = -2
call DMUMPS(mumps_par)

if (mumps_par%INFOG(1).lt.0) then
    write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",                           &
                            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
    goto 500
endif

500 return
!--------------------------------------------------------------------------! 
end subroutine mumps_sub
