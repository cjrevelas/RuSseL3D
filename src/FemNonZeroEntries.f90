!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine FemNonZeroEntries(stepSize, nodeBelongsToDirichletFace, elemcon)
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod,         only: F_m, A_m, numNonZeroEntries
use geometry_mod,    only: numTotalNodePairs, numNodesLocalTypeDomain, numElementsTypeDomain, &
                           numNodes, nodePairingXXhash, nodePairingYYhash, nodePairingZZhash
use constants_mod,   only: tol
use parser_vars_mod, only: periodicAxisId, periodicity, mumpsMatrixType
use flags_mod,       only: mumpsAsymmetric, mumpsPositiveDefinite, mumpsGeneralSymmetric
use fhash_module__ints_double
!#define PRINT_AFULL
#ifdef PRINT_AFULL
use iofiles_mod, only: IO_A_matrix_full
#endif
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
logical, intent(in), dimension(numNodes)                          :: nodeBelongsToDirichletFace
logical, dimension(numNodesLocalTypeDomain*numElementsTypeDomain) :: setDiagToOne

real(8), intent(in) :: stepSize

type(fhash_type__ints_double), intent(inout) :: elemcon

integer :: nodePair, row, col

#ifdef PRINT_AFULL
integer :: entry

real(8), allocatable, dimension(:,:) :: A_full
#endif
!------------------------------------------------------------------------------------------------------!
F_m%g  = F_m%c + stepSize * (F_m%k + F_m%w)
F_m%rh = F_m%c

! Apply periodic boundary conditions
if (periodicAxisId(1)) call FemApplyPeriodicity(nodePairingXXhash, elemcon)
if (periodicAxisId(2)) call FemApplyPeriodicity(nodePairingYYhash, elemcon)
if (periodicAxisId(3)) call FemApplyPeriodicity(nodePairingZZhash, elemcon)

if (periodicity.eq.2) then
  if (periodicAxisId(1).AND.periodicAxisId(2)) call FemPeriodicEdges(elemcon)
  if (periodicAxisId(1).AND.periodicAxisId(3)) call FemPeriodicEdges(elemcon)
  if (periodicAxisId(2).AND.periodicAxisId(3)) call FemPeriodicEdges(elemcon)
endif

if (periodicity.eq.3) then
  call FemPeriodicEdges(elemcon)
  call FemPeriodicCorners(elemcon)
endif

! Prepare stiffness matrix for Dirichlet boundary conditions
! In case the matrix is symmetric, remove the zero lines and rows diagonal componets with Dirichlet BC q=0.
setDiagToOne = .True.
if ((mumpsMatrixType.eq.mumpsPositiveDefinite).or.(mumpsMatrixType.eq.mumpsGeneralSymmetric)) then
  do nodePair = 1, numTotalNodePairs
    if (F_m%isZero(nodePair)) cycle
    if (F_m%row(nodePair)==0) cycle

    row = F_m%row(nodePair)
    col = F_m%col(nodePair)

    if (row > col) F_m%g(nodePair) = 0.0d0

    if (nodeBelongsToDirichletFace(row).OR.nodeBelongsToDirichletFace(col)) then
      F_m%g(nodePair) = 0.0d0
      if ((row==col).AND.setDiagToOne(row)) then
        F_m%g(nodePair)   = 1.0d0
        setDiagToOne(row) = .false.
      endif
    endif
  enddo
endif

if (mumpsMatrixType.eq.mumpsAsymmetric) then
  do nodePair = 1, numTotalNodePairs
    if (F_m%isZero(nodePair)) cycle

    if (F_m%row(nodePair)==0) cycle

    row = F_m%row(nodePair)
    col = F_m%col(nodePair)

    if (nodeBelongsToDirichletFace(row)) then
      F_m%g(nodePair) = 0.0d0
      if ((row==col).AND.setDiagToOne(row)) then
        F_m%g(nodePair)   =  1.0d0
        setDiagToOne(row) = .False.
      endif
    endif
  enddo
endif

! Determine non_zero entries
numNonZeroEntries = 0
do nodePair = 1, numTotalNodePairs
  if (ABS(F_m%g(nodePair)) > tol) numNonZeroEntries = numNonZeroEntries + 1
enddo

allocate(A_m%value(numNonZeroEntries))
allocate(A_m%col(numNonZeroEntries))
allocate(A_m%row(numNonZeroEntries))

numNonZeroEntries = 0
do nodePair = 1, numTotalNodePairs
  if (ABS(F_m%g(nodePair)) > tol) then
    numNonZeroEntries = numNonZeroEntries + 1

    A_m%value(numNonZeroEntries) = F_m%g(nodePair)
    A_m%row(numNonZeroEntries)   = F_m%row(nodePair)
    A_m%col(numNonZeroEntries)   = F_m%col(nodePair)
  endif
enddo

#ifdef PRINT_AFULL
allocate(A_full(numNodes,numNodes))
A_full = 0.0d0

do entry = 1, numNonZeroEntries
  row = A_m%row(entry)
  col = A_m%col(entry)

  A_full(row,col) = A_m%value(entry)
enddo

open(unit=255, file = IO_A_matrix_full)
do row = 1, numNodes
  write(255,*) (A_full(row,col), col = 1, numNodes)
enddo
close(255)

deallocate(A_full)
#endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine FemNonZeroEntries
