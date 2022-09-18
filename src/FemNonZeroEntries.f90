!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine FemNonZeroEntries(ds, mumpsMatrixType, nodeBelongsToDirichletFace, elemcon)
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod,         only: F_m, A_m, NNZ
use geometry_mod,    only: numTotalNodePairs, numNodesLocalTypeDomain, numElementsTypeDomain, &
                           numNodes, nodePairingXXhash, nodePairingYYhash, nodePairingZZhash, &
                           nodePairingYYhashInverse, nodePairingZZhashInverse
use constants_mod,   only: tol
use parser_vars_mod, only: periodicAxisId, periodicity
use flags_mod,       only: mumps_asymm, mumps_posDef, mumps_genSymm
use fhash_module__ints_double
!#define PRINT_AFULL
#ifdef PRINT_AFULL
use iofiles_mod, only: A_matrix_full
#endif
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: mumpsMatrixType
integer             :: ii, jj, kk

logical, intent(in), dimension(numNodes)                          :: nodeBelongsToDirichletFace
logical, dimension(numNodesLocalTypeDomain*numElementsTypeDomain) :: set_diag_to_one

real(8), intent(in) :: ds

type(fhash_type__ints_double), intent(inout) :: elemcon

#ifdef PRINT_AFULL
real(8), allocatable, dimension(:,:) :: A_full
#endif
!------------------------------------------------------------------------------------------------------!
F_m%g  = F_m%c + ds * (F_m%k + F_m%w)
F_m%rh = F_m%c

! Apply periodic boundary conditions
if (periodicAxisId(1)) call FemApplyPeriodicity(nodePairingXXhash, elemcon)
if (periodicAxisId(2)) call FemApplyPeriodicity(nodePairingYYhash, elemcon)
if (periodicAxisId(3)) call FemApplyPeriodicity(nodePairingZZhash, elemcon)

if (periodicity.eq.2) then
  if (periodicAxisId(1).AND.periodicAxisId(2)) call FemPeriodicEdges(nodePairingXXhash, nodePairingYYhash, nodePairingYYhashInverse, elemcon)
  if (periodicAxisId(1).AND.periodicAxisId(3)) call FemPeriodicEdges(nodePairingXXhash, nodePairingZZhash, nodePairingZZhashInverse, elemcon)
  if (periodicAxisId(2).AND.periodicAxisId(3)) call FemPeriodicEdges(nodePairingYYhash, nodePairingZZhash, nodePairingZZhashInverse, elemcon)
endif

if (periodicity.eq.3) then
  call FemPeriodicEdges(nodePairingXXhash, nodePairingYYhash, nodePairingYYhashInverse, elemcon)
  call FemPeriodicCorners(elemcon)
endif

! Prepare stiffness matrix for Dirichlet boundary conditions
! In case the matrix is symmetric, remove the zero lines and rows diagonal componets with Dirichlet BC q=0.
set_diag_to_one=.true.
if ((mumpsMatrixType.eq.mumps_posDef).or.(mumpsMatrixType.eq.mumps_genSymm)) then
  do kk = 1, numTotalNodePairs
    if (F_m%is_zero(kk)) cycle
    if (F_m%row(kk)==0)  cycle

    ii = F_m%row(kk)
    jj = F_m%col(kk)

    if (ii > jj) F_m%g(kk) = 0.0d0

    if (nodeBelongsToDirichletFace(ii).OR.nodeBelongsToDirichletFace(jj)) then
      F_m%g(kk) = 0.0d0
      if (ii==jj.and.set_diag_to_one(ii)) then
        F_m%g(kk)           = 1.0d0
        set_diag_to_one(ii) = .false.
      endif
    endif
  enddo
endif

if (mumpsMatrixType.eq.mumps_asymm) then
  do kk = 1, numTotalNodePairs
    if (F_m%is_zero(kk)) cycle

    if (F_m%row(kk)==0) cycle

    jj = F_m%col(kk)
    ii = F_m%row(kk)

    if (nodeBelongsToDirichletFace(ii)) then
      F_m%g(kk) = 0.0d0
      if (ii==jj.and.set_diag_to_one(ii)) then
        F_m%g(kk)           =  1.0d0
        set_diag_to_one(ii) = .false.
      endif
    endif
  enddo
endif

! Determine non_zero entries
NNZ = 0
do kk = 1, numTotalNodePairs
  if (ABS(F_m%g(kk)) > tol) NNZ = NNZ + 1
enddo

allocate(A_m%value(NNZ))
allocate(A_m%col(NNZ))
allocate(A_m%row(NNZ))

NNZ = 0
do kk = 1, numTotalNodePairs
  if (ABS(F_m%g(kk)) > tol) then
    NNZ = NNZ + 1

    A_m%value(NNZ) = F_m%g(kk)
    A_m%row(NNZ)   = F_m%row(kk)
    A_m%col(NNZ)   = F_m%col(kk)
  endif
enddo

#ifdef PRINT_AFULL
allocate(A_full(numNodes,numNodes))
A_full = 0.0d0

do kk = 1, NNZ
  ii = A_m%row(kk)
  jj = A_m%col(kk)

  A_full(ii,jj) = A_m%value(kk)
enddo

open(unit=255, file = A_matrix_full)
do ii = 1, numNodes
  write(255,*) (A_full(ii,jj), jj = 1, numNodes)
enddo
close(255)

deallocate(A_full)
#endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine FemNonZeroEntries
