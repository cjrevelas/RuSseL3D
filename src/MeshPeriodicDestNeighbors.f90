subroutine MeshPeriodicDestNeighbors(isDestPeriodicNodeXX, isDestPeriodicNodeYY, isDestPeriodicNodeZZ, elemcon, nodePairId)
!----------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use geometry_mod,    only: numNodes, numBulkNodePairs, numTotalNodePairs, &
                           nodePairingXXhashInverse,                      &
                           nodePairingYYhashInverse,                      &
                           nodePairingZZhashInverse
use parser_vars_mod, only: periodicAxisId
use kcw_mod,         only: F_m
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: elemcon

logical, dimension(numNodes), intent(in) :: isDestPeriodicNodeXX
logical, dimension(numNodes), intent(in) :: isDestPeriodicNodeYY
logical, dimension(numNodes), intent(in) :: isDestPeriodicNodeZZ

integer, dimension(numTotalNodePairs), intent(inout) :: nodePairId

type(ints_type) :: elemconKey, nodePairingXXkeyInverse, nodePairingYYkeyInverse, nodePairingZZkeyInverse

integer :: node1, node2, pair, pairId
integer :: nodePairingXXvalueInverse, nodePairingYYvalueInverse, nodePairingZZvalueInverse
integer :: elemconValueBulk

logical :: success
!----------------------------------------------------------------------------------------------------------------------------!
allocate(elemconKey%ints(2))

if (periodicAxisId(1)) allocate(nodePairingXXkeyInverse%ints(1))
if (periodicAxisId(2)) allocate(nodePairingYYkeyInverse%ints(1))
if (periodicAxisId(3)) allocate(nodePairingZZkeyInverse%ints(1))

do pair = 1, numBulkNodePairs
  node1  = F_m%row(pair)     ! 24  -> neumann
  node2  = F_m%col(pair)     ! 23  -> neumann
  pairId = nodePairId(pair)  ! 103 = ii = nodePairId(ii)

  if ((node1==node2).and.(isDestPeriodicNodeXX(node1)).and.(isDestPeriodicNodeXX(node2))) cycle
  if ((node1==node2).and.(isDestPeriodicNodeYY(node1)).and.(isDestPeriodicNodeYY(node2))) cycle
  if ((node1==node2).and.(isDestPeriodicNodeZZ(node1)).and.(isDestPeriodicNodeZZ(node2))) cycle

  if (isDestPeriodicNodeXX(node1)) then
    nodePairingXXkeyInverse%ints(1) = node1
    call nodePairingXXhashInverse%get(nodePairingXXkeyInverse, nodePairingXXvalueInverse)
    F_m%row(pair) = nodePairingXXvalueInverse
    node1         = nodePairingXXvalueInverse ! This change is necessary, so that we know if a src along the X dimension is a dst along another dimension
  endif

  if (isDestPeriodicNodeXX(node2)) then
    nodePairingXXkeyInverse%ints(1) = node2
    call nodePairingXXhashInverse%get(nodePairingXXkeyInverse, nodePairingXXvalueInverse)
    F_m%col(pair) = nodePairingXXvalueInverse
    node2         = nodePairingXXvalueInverse
  endif

  if (isDestPeriodicNodeYY(node1)) then
    nodePairingYYkeyInverse%ints(1) = node1
    call nodePairingYYhashInverse%get(nodePairingYYkeyInverse, nodePairingYYvalueInverse)
    F_m%row(pair) = nodePairingYYvalueInverse
    node1         = nodePairingYYvalueInverse
  endif

  if (isDestPeriodicNodeYY(node2)) then
    nodePairingYYkeyInverse%ints(1) = node2
    call nodePairingYYhashInverse%get(nodePairingYYkeyInverse, nodePairingYYvalueInverse)
    F_m%col(pair) = nodePairingYYvalueInverse
    node2         = nodePairingYYvalueInverse
  endif

  if (isDestPeriodicNodeZZ(node1)) then
    nodePairingZZkeyInverse%ints(1) = node1
    call nodePairingZZhashInverse%get(nodePairingZZkeyInverse, nodePairingZZvalueInverse)
    F_m%row(pair) = nodePairingZZvalueInverse
    node1         = nodePairingZZvalueInverse
  endif

  if (isDestPeriodicNodeZZ(node2)) then
    nodePairingZZkeyInverse%ints(1) = node2
    call nodePairingZZhashInverse%get(nodePairingZZkeyInverse, nodePairingZZvalueInverse)
    F_m%col(pair) = nodePairingZZvalueInverse
    node2         = nodePairingZZvalueInverse
  endif

  elemconKey%ints(1) = F_m%row(pair) ! 24 -> periodic
  elemconKey%ints(2) = F_m%col(pair) ! 7  -> periodic

  call elemcon%get(elemconKey, elemconValueBulk, success) ! hashBulk(24,7) = 201

  if (success) then
    if (elemconValueBulk > pairId) then
      call elemcon%set(elemconKey, pairId)  ! elemconValueBulk will become equal to pairId for next instances of the same pair,   -> hash(24,7)      = 103
      nodePairId(pair) = pairId             ! nodePairId will become equal to pairId (if is not already)                          -> nodePairId(103) = 103
    else
      !call elemcon%set(elemconKey, elemconValueBulk) ! THIS LINE MAYBE IS NOT NEEDED
      nodePairId(pair) = elemconValueBulk             ! nodePairId(201) = 103
    endif
  else
    call elemcon%set(elemconKey, pairId) ! hash(24,7)     = 103
    nodePairId(pair) = pairId            ! nodePairId(ii) = 103
  endif
enddo

deallocate(elemconKey%ints)

if (periodicAxisId(1)) deallocate(nodePairingXXkeyInverse%ints)
if (periodicAxisId(2)) deallocate(nodePairingYYkeyInverse%ints)
if (periodicAxisId(3)) deallocate(nodePairingZZkeyInverse%ints)

return
!----------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicDestNeighbors
