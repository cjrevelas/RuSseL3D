!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine mesh_build_node_pairing(globalNodeIdTypeFace, axis, face1_hash, face2_hash, node_pairing_hash)
!----------------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module

use geometry_mod, only: nodeCoord, numNodesLocalTypeFace, numElementsTypeFace
use iofiles_mod,  only: IO_nodePairingXX, IO_nodePairingYY, IO_nodePairingZZ
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
character(len=1), intent(in) :: axis

integer, dimension(numNodesLocalTypeFace, numElementsTypeFace), intent(in) :: globalNodeIdTypeFace

type(fhash_type__ints_double), intent(out) :: node_pairing_hash
type(fhash_type_iterator__ints_double)     :: node_pairing_it
type(ints_type)                            :: node_pairing_key
integer                                    :: node_pairing_value

type(fhash_type__ints_double), intent(inout) :: face1_hash, face2_hash
type(fhash_type_iterator__ints_double)       :: face1_it, face2_it
type(ints_type)                              :: face1_key, face2_key
integer                                      :: face1_value, face2_value

integer :: ii, jj, kk, mm
integer :: elem1, elem2, node1, node2
integer :: dimension1=0, dimension2=0
!----------------------------------------------------------------------------------------------------------------------------------!
if (axis=='x') then
  dimension1 = 2 ! y
  dimension2 = 3 ! z

  open(unit=1111, file=IO_nodePairingXX)
elseif (axis=='y') then
  dimension1 = 1 ! x
  dimension2 = 3 ! z

  open(unit=1111, file=IO_nodePairingYY)
elseif (axis=='z') then
  dimension1 = 1 ! x
  dimension2 = 2 !y

  open(unit=1111, file=IO_nodePairingZZ)
endif

allocate(node_pairing_key%ints(1))
call node_pairing_hash%reserve(face1_hash%key_count())
call face1_it%begin(face1_hash)
do kk = 1, face1_hash%key_count()
  call face1_it%next(face1_key, face1_value)
  elem1 = face1_key%ints(1)

  call face2_it%begin(face2_hash)
  do mm = 1, face2_hash%key_count()
    call face2_it%next(face2_key, face2_value)
    elem2 = face2_key%ints(1)
    do ii = 1, 3
      node1 = globalNodeIdTypeFace(ii,elem1)
      do jj = 1, 3
        node2 = globalNodeIdTypeFace(jj, elem2)
        if ((ABS(nodeCoord(dimension1,node1)-nodeCoord(dimension1,node2))<1.0d-12).AND.(ABS(nodeCoord(dimension2,node1)-nodeCoord(dimension2,node2))<1.0d-12)) then
          node_pairing_key%ints(1) = node1
          call node_pairing_hash%set(node_pairing_key, node2)
        endif
      enddo
    enddo
  enddo
enddo

#ifdef DEBUG_OUTPUTS
call node_pairing_it%begin(node_pairing_hash)

do kk = 1, node_pairing_hash%key_count()
  call node_pairing_it%next(node_pairing_key, node_pairing_value)
  write(1111,*) node_pairing_key%ints(1), node_pairing_value
enddo

close(1111)
#endif

deallocate(node_pairing_key%ints)
deallocate(face1_key%ints)
deallocate(face2_key%ints)

call face1_hash%clear()
call face2_hash%clear()

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine mesh_build_node_pairing
