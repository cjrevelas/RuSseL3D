!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshPeriodicNodePairs(globalNodeIdTypeFace, axis, faceOneHash, faceTwoHash, isDestPeriodicNode, nodePairingHash, nodePairingHashInverse)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module

use geometry_mod,  only: nodeCoord, numNodesLocalTypeFace, numElementsTypeFace, numNodes
use iofiles_mod,   only: IO_nodePairingXX, IO_nodePairingYY, IO_nodePairingZZ
use constants_mod, only: tol
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
character(len=1), intent(in) :: axis

integer, intent(in), dimension(numNodesLocalTypeFace, numElementsTypeFace) :: globalNodeIdTypeFace

logical, intent(out), dimension(numNodes) :: isDestPeriodicNode

type(fhash_type__ints_double), intent(out) :: nodePairingHash, nodePairingHashInverse
type(fhash_type_iterator__ints_double)     :: nodePairingIt, nodePairingItInverse
type(ints_type)                            :: nodePairingKey, nodePairingKeyInverse
integer                                    :: nodePairingValue, nodePairingValueInverse

type(fhash_type__ints_double), intent(inout) :: faceOneHash, faceTwoHash
type(fhash_type_iterator__ints_double)       :: faceOneIt, faceTwoIt
type(ints_type)                              :: faceOneKey, faceTwoKey
integer                                      :: faceOneValue, faceTwoValue

integer :: localNodeIndex1, localNodeIndex2
integer :: globalNodeId1, globalNodeId2
integer :: keyIndex1, keyIndex2
integer :: element1, element2
integer :: dimension1 = 0, dimension2 = 0
!----------------------------------------------------------------------------------------------------------------------------------!
if (axis=='x') then
  dimension1 = 2 ! y
  dimension2 = 3 ! z

#ifdef DEBUG_OUTPUTS
  open(unit=1111, file=IO_nodePairingXX)
#endif
elseif (axis=='y') then
  dimension1 = 1 ! x
  dimension2 = 3 ! z

#ifdef DEBUG_OUTPUTS
  open(unit=1111, file=IO_nodePairingYY)
#endif
elseif (axis=='z') then
  dimension1 = 1 ! x
  dimension2 = 2 ! y

#ifdef DEBUG_OUTPUTS
  open(unit=1111, file=IO_nodePairingZZ)
#endif
endif

allocate(nodePairingKey%ints(1))
allocate(nodePairingKeyInverse%ints(1))

CALL nodePairingHash%reserve(faceOneHash%key_count())
CALL nodePairingHashInverse%reserve(faceTwoHash%key_count())

isDestPeriodicNode = .False.

CALL faceOneIt%begin(faceOneHash)
do keyIndex1 = 1, faceOneHash%key_count()
  CALL faceOneIt%next(faceOneKey, faceOneValue)
  element1 = faceOneKey%ints(1)

  CALL faceTwoIt%begin(faceTwoHash)
  do keyIndex2 = 1, faceTwoHash%key_count()
    CALL faceTwoIt%next(faceTwoKey, faceTwoValue)
    element2 = faceTwoKey%ints(1)
    do localNodeIndex1 = 1, 3
      globalNodeId1 = globalNodeIdTypeFace(localNodeIndex1,element1)
      do localNodeIndex2 = 1, 3
        globalNodeId2 = globalNodeIdTypeFace(localNodeIndex2, element2)
        if ((ABS(nodeCoord(dimension1,globalNodeId1) - nodeCoord(dimension1,globalNodeId2)) < tol).AND.&
            (ABS(nodeCoord(dimension2,globalNodeId1) - nodeCoord(dimension2,globalNodeId2)) < tol)) then
          nodePairingKey%ints(1) = globalNodeId1
          CALL nodePairingHash%set(nodePairingKey, globalNodeId2)

          nodePairingKeyInverse%ints(1) = globalNodeId2
          CALL nodePairingHashInverse%set(nodePairingKeyInverse, globalNodeId1)

          isDestPeriodicNode(globalNodeId2) = .True.
        endif
      enddo
    enddo
  enddo
enddo

#ifdef DEBUG_OUTPUTS
CALL nodePairingIt%begin(nodePairingHash)
CALL nodePairingItInverse%begin(nodePairingHashInverse)

do keyIndex1 = 1, nodePairingHash%key_count()
  CALL nodePairingIt%next(nodePairingKey, nodePairingValue)
  write(1111,*) nodePairingKey%ints(1), nodePairingValue
enddo

write(1111,*)

do keyIndex1 = 1, nodePairingHashInverse%key_count()
  CALL nodePairingItInverse%next(nodePairingKeyInverse, nodePairingValueInverse)
  write(1111,*) nodePairingKeyInverse%ints(1), nodePairingValueInverse
enddo

close(1111)
#endif

deallocate(nodePairingKey%ints)
deallocate(nodePairingKeyInverse%ints)
deallocate(faceOneKey%ints)
deallocate(faceTwoKey%ints)

CALL faceOneHash%clear()
CALL faceTwoHash%clear()

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicNodePairs
