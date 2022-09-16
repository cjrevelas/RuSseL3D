!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshPeriodicNodePairs(globalNodeIdTypeFace, axis, faceOneHash, faceTwoHash, isDestPeriodicNode, nodePairingHash, nodePairingHashInverse)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module

use geometry_mod, only: nodeCoord, numNodesLocalTypeFace, numElementsTypeFace, numNodes
use iofiles_mod,  only: IO_nodePairingXX, IO_nodePairingYY, IO_nodePairingZZ
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

integer :: ii, jj, kk, mm
integer :: elem1, elem2, node1, node2
integer :: dimension1=0, dimension2=0
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

call nodePairingHash%reserve(faceOneHash%key_count())
call nodePairingHashInverse%reserve(faceTwoHash%key_count())

isDestPeriodicNode = .false.

call faceOneIt%begin(faceOneHash)
do kk = 1, faceOneHash%key_count()
  call faceOneIt%next(faceOneKey, faceOneValue)
  elem1 = faceOneKey%ints(1)

  call faceTwoIt%begin(faceTwoHash)
  do mm = 1, faceTwoHash%key_count()
    call faceTwoIt%next(faceTwoKey, faceTwoValue)
    elem2 = faceTwoKey%ints(1)
    do ii = 1, 3
      node1 = globalNodeIdTypeFace(ii,elem1)
      do jj = 1, 3
        node2 = globalNodeIdTypeFace(jj, elem2)
        if ((ABS(nodeCoord(dimension1,node1)-nodeCoord(dimension1,node2))<1.0d-12).AND.(ABS(nodeCoord(dimension2,node1)-nodeCoord(dimension2,node2))<1.0d-12)) then
          nodePairingKey%ints(1) = node1
          call nodePairingHash%set(nodePairingKey, node2)

          nodePairingKeyInverse%ints(1) = node2
          call nodePairingHashInverse%set(nodePairingKeyInverse, node1)

          isDestPeriodicNode(node2) = .True.
        endif
      enddo
    enddo
  enddo
enddo

#ifdef DEBUG_OUTPUTS
call nodePairingIt%begin(nodePairingHash)
call nodePairingItInverse%begin(nodePairingHashInverse)

do kk = 1, nodePairingHash%key_count()
  call nodePairingIt%next(nodePairingKey, nodePairingValue)
  write(1111,*) nodePairingKey%ints(1), nodePairingValue
enddo

write(1111,*)

do kk = 1, nodePairingHashInverse%key_count()
  call nodePairingItInverse%next(nodePairingKeyInverse, nodePairingValueInverse)
  write(1111,*) nodePairingKeyInverse%ints(1), nodePairingValueInverse
enddo

close(1111)
#endif

deallocate(nodePairingKey%ints)
deallocate(nodePairingKeyInverse%ints)
deallocate(faceOneKey%ints)
deallocate(faceTwoKey%ints)

call faceOneHash%clear()
call faceTwoHash%clear()

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicNodePairs
