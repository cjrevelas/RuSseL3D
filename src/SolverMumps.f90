!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine SolverMumps(mumpsMatrixType)
!--------------------------------------------------------------------------!
use kcw_mod,      only: A_m, rdiag1, numNonZeroEntries
use geometry_mod, only: numNodes
use flags_mod,    only: mumpsAsymmetric, mumpsPositiveDefinite, &
                        mumpsGeneralSymmetric
use error_handling_mod
#ifdef USEMPI
use mpistuff_mod
#endif
!--------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------!
include "mpif.h"
include "dmumps_struc.h"
!--------------------------------------------------------------------------!
type(DMUMPS_STRUC) :: mumps_par

integer, intent(in) :: mumpsMatrixType
integer             :: node
!--------------------------------------------------------------------------!
! Define a communicator for the package
mumps_par%COMM = MPI_COMM_WORLD

! Initialize an instance of the package for LU-factorization
mumps_par%PAR = 1  ! Working host processor

! Set the type of the matrix
if (mumpsMatrixType.eq.mumpsAsymmetric) then
  mumps_par%SYM = 0
elseif (mumpsMatrixType.eq.mumpsPositiveDefinite) then
  mumps_par%SYM       = 1
  mumps_par%ICNTL(13) = 0
elseif (mumpsMatrixType.eq.mumpsGeneralSymmetric) then
  mumps_par%SYM     = 2
  mumps_par%CNTL(1) = 0
else
  ERROR_MESSAGE="MUMPS SUBROUTINE: mumpsMatrixType not between 0-2."
  CALL exitWithError(1, 2, 1, ERROR_MESSAGE)
endif

mumps_par%JOB = -1

CALL DMUMPS(mumps_par)

if (mumps_par%INFOG(1).lt.0) then
  write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                                              "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
  goto 500
endif

#ifndef MUMPS_REPORT
mumps_par%ICNTL(1) = -1
mumps_par%ICNTL(2) = -1
mumps_par%ICNTL(3) = -1
mumps_par%ICNTL(4) = -1
#endif

! Set mumps options here
!mumps_par%ICNTL(28) = 2 !Parallel Ordering tools
!mumps_par%ICNTL(7)  = 4 !Choose ordering scheme

! Define problem on the master process (id = 0)
if (mumps_par%MYID.eq.0) then
  mumps_par%N   = numNodes
  mumps_par%NNZ = numNonZeroEntries

  allocate(mumps_par%IRN(mumps_par%NNZ))
  allocate(mumps_par%JCN(mumps_par%NNZ))
  allocate(mumps_par%A(mumps_par%NNZ))
  allocate(mumps_par%RHS(mumps_par%N))

  mumps_par%IRN = A_m%row
  mumps_par%JCN = A_m%col
  mumps_par%A   = A_m%value

  do node = 1, mumps_par%N
    mumps_par%RHS(node) = rdiag1(node)
  enddo
endif

! Call package for matrix factorization and solution
mumps_par%JOB = 6

CALL DMUMPS(mumps_par)

if (mumps_par%INFOG(1).lt.0) then
  write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                                              "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
  goto 500
endif

! Solution has been assembled on the master process
if (mumps_par%MYID.eq.0) then
  do node = 1, mumps_par%N
    rdiag1(node) = mumps_par%RHS(node)
  enddo
endif

! Deallocate user data
if (mumps_par%MYID.eq.0) then
  deallocate(mumps_par%IRN)
  deallocate(mumps_par%JCN)
  deallocate(mumps_par%A)
  deallocate(mumps_par%RHS)
endif

! Destruct the instance (deallocate internal data structures)
mumps_par%JOB = -2
CALL DMUMPS(mumps_par)

if (mumps_par%INFOG(1).lt.0) then
  write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                                              "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
  goto 500
endif

500 return
!--------------------------------------------------------------------------!
end subroutine SolverMumps
