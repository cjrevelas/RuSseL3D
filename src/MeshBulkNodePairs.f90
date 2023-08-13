!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshBulkNodePairs(elemcon)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use geometry_mod, only: numNodesLocalTypeDomain, numElementsTypeDomain, &
                           numTotalNodePairs, nodePairId, globalNodeIdTypeDomain
use kcw_mod,      only: F_m
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
logical :: success

type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemconKey
integer                                      :: elemconNumOfKeys, elemconValue

integer :: localNodeIndex1, localNodeIndex2, element
integer :: nodePair
!----------------------------------------------------------------------------------------------------------------------------------!
allocate(nodePairId(numTotalNodePairs))
nodePairId = 0

! Assembly the nodePairId (hash) matrix
allocate(elemconKey%ints(2)) ! Each key is defined by a pair (2) of nodes

! Total number of required keys
elemconNumOfKeys = 2 * numNodesLocalTypeDomain * numElementsTypeDomain
call elemcon%reserve(elemconNumOfKeys)

nodePair = 0
do element = 1, numElementsTypeDomain
  do localNodeIndex1 = 1, numNodesLocalTypeDomain
    do localNodeIndex2 = 1, numNodesLocalTypeDomain
      nodePair = nodePair + 1

      ! Define the pair of nodes to be examined and assigned a elemconValue
      elemconKey%ints(1) = globalNodeIdTypeDomain(localNodeIndex1,element)
      elemconKey%ints(2) = globalNodeIdTypeDomain(localNodeIndex2,element)

      F_m%row(nodePair) = elemconKey%ints(1)
      F_m%col(nodePair) = elemconKey%ints(2)

      ! Assign elemconValue to the pair
      call elemcon%get(elemconKey, elemconValue, success)

      if (success) then
        nodePairId(nodePair) = elemconValue  ! This pair has already been met, thus assigned an elemconValue
      else
        call elemcon%set(elemconKey, nodePair) ! Store the new elemconValue for next iteration's check
        nodePairId(nodePair) = nodePair       ! This pair is met for the first time
      endif
    enddo
  enddo
enddo

deallocate(elemconKey%ints)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshBulkNodePairs
