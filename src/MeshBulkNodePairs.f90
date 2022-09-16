!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshBulkNodePairs(elemcon)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use geometry_mod,    only: numNodesLocalTypeDomain, numElementsTypeDomain, numTotalNodePairs, nodePairId, globalNodeIdTypeDomain
use kcw_mod,         only: F_m
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
logical :: success

type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemconKey
integer                                      :: elemconNumOfKeys, elemconValue

integer :: ii, jj, mm
integer :: node_pair
!----------------------------------------------------------------------------------------------------------------------------------!
allocate(nodePairId(numTotalNodePairs))
nodePairId = 0

! Assembly the nodePairId (hash) matrix
allocate(elemconKey%ints(2)) ! Each key is defined by a pair (2) of nodes

! Total number of required keys
elemconNumOfKeys = 2 * numNodesLocalTypeDomain * numElementsTypeDomain
call elemcon%reserve(elemconNumOfKeys)

node_pair = 0
do mm = 1, numElementsTypeDomain
  do jj = 1, numNodesLocalTypeDomain
    do ii = 1, numNodesLocalTypeDomain
      node_pair = node_pair + 1

      ! Define the pair of nodes to be examined and assigned a elemconValue
      elemconKey%ints(1) = globalNodeIdTypeDomain(jj,mm)
      elemconKey%ints(2) = globalNodeIdTypeDomain(ii,mm)

      F_m%row(node_pair) = elemconKey%ints(1)
      F_m%col(node_pair) = elemconKey%ints(2)

      ! Assign elemconValue to the pair
      call elemcon%get(elemconKey, elemconValue, success)

      if (success) then
        nodePairId(node_pair) = elemconValue  ! This pair has already been met, thus assigned a elemconValue
      else
        call elemcon%set(elemconKey, node_pair) ! Store the new elemconValue for next iteration's check
        nodePairId(node_pair) = node_pair        ! This pair is met for the first time
      endif
    enddo
  enddo
enddo

deallocate(elemconKey%ints)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshBulkNodePairs
