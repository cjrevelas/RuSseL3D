
!This file is part of MUMPS 5.2.1, released
!on Fri Jun 14 14:46:05 UTC2019

      SUBROUTINE MUMPS_SUB(numnp,rdiag,u1,i9)
      use kcw
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC) mumps_par
      INTEGER EX_INPUT
        real*8 rdiag(1:1000)
      INTEGER IERR, I, numnp
      INTEGER(8) I8,i9
      double precision ::u1(numnp)


      

      CALL MPI_INIT(IERR)
!Define a communicator for the package.
      mumps_par%COMM = MPI_COMM_WORLD
!Initialize an instance of the package
!for L U factorization (sym = 0, with working host)
      mumps_par%JOB = -1
      mumps_par%SYM = 0
      mumps_par%PAR = 1
      CALL DMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
                 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
       GOTO 500
      END IF
      !tolis
       do i=1,4
      mumps_par%ICNTL(i)=-1
      end do
      !tolis
!Define problem on the host (processor 0)
      IF ( mumps_par%MYID .eq. 0 ) THEN
!        READ(5,*) mumps_par%N


        mumps_par%N=numnp
        mumps_par%NNZ=non_zero

!        READ(5,*) mumps_par%NNZ
        ALLOCATE( mumps_par%IRN ( mumps_par%NNZ ) )
        ALLOCATE( mumps_par%JCN ( mumps_par%NNZ ) )
        ALLOCATE( mumps_par%A( mumps_par%NNZ ) )
        ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
 !       allocate(u1(mumps_par%N))
        mumps_par%IRN=c_m%row
        mumps_par%JCN=c_m%col
        mumps_par%A=c_m%value
        
         
         
        DO I = 1, mumps_par%N
           mumps_par%RHS(I)=rdiag(i)
        END DO
      END IF
      
      
!DEBUG       
        if(i9==0)then
         open(file='mumpsinside.txt',unit=55)
             DO I8 = 1, mumps_par%NNZ
                 write(55,*) mumps_par%IRN(I8),mumps_par%JCN(I8),mumps_par%A(I8)
              END DO
              
              DO I8 = 1, mumps_par%N
                 write(55,*) mumps_par%RHS(I8)
              END DO
              
         end if 
!DEBUG

!Call package for solution
      mumps_par%JOB = 6
      CALL DMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
                 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
       GOTO 500
      END IF
!Solution has been assembled on the host
        if ( mumps_par%MYID .eq. 0 )then
        open(file='u1.txt',unit=53)
            do i=1,mumps_par%N
              u1(i) =mumps_par%RHS(i)
                 write(53,*)i,u1(i)
            end do
        end if
!Deallocate user data


      IF ( mumps_par%MYID .eq. 0 )THEN
        DEALLOCATE( mumps_par%IRN )
        DEALLOCATE( mumps_par%JCN )
        DEALLOCATE( mumps_par%A   )
        DEALLOCATE( mumps_par%RHS )
      END IF
!Destroy the instance (deallocate internal data structures)
      mumps_par%JOB = -2
      CALL DMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
                 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
       GOTO 500
      END IF
 500  CALL MPI_FINALIZE(IERR)
      END SUBROUTINE
