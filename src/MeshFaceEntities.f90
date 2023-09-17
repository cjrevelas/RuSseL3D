!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshFaceEntities(axis, boundaryCoord, globalNodeIdTypeFace, faceEntityHash, faceHash)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use geometry_mod,  only: numNodesLocalTypeFace, numElementsTypeFace, nodeCoord
use constants_mod, only: tol
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
integer, intent(in) :: axis
real(8), intent(in) :: boundaryCoord

integer, intent(in), dimension(numNodesLocalTypeFace, numElementsTypeFace) :: globalNodeIdTypeFace

type(fhash_type__ints_double), intent(inout) :: faceEntityHash
type(fhash_type_iterator__ints_double)       :: faceEntityIt
type(ints_type)                              :: faceEntityKey
integer                                      :: faceEntityValue

type(fhash_type__ints_double), intent(out) :: faceHash
type(ints_type)                            :: faceKey
integer                                    :: faceSize

real(8) :: centerOfMassCoord

integer :: keyIndex, localNodeId, globalNodeId, element
!----------------------------------------------------------------------------------------------------------------------------------!
faceSize = 0

! Determine the size of face entity hashes for proper allocation
CALL faceEntityIt%begin(faceEntityHash)
do keyIndex = 1, faceEntityHash%key_count()
  centerOfMassCoord = 0.d0

  CALL faceEntityIt%next(faceEntityKey, faceEntityValue)

  element = faceEntityKey%ints(1)

  do localNodeId = 1, 3
    globalNodeId = globalNodeIdTypeFace(localNodeId,element)

    centerOfMassCoord = centerOfMassCoord + nodeCoord(axis,globalNodeId)
  enddo

  centerOfMassCoord = centerOfMassCoord / 3.0d0

  if (ABS(centerOfMassCoord - boundaryCoord) < tol) faceSize = faceSize + 1
enddo

allocate(faceKey%ints(1))

CALL faceHash%reserve(faceSize)

! Fill the face entity hashes
CALL faceEntityIt%begin(faceEntityHash)
do keyIndex = 1, faceEntityHash%key_count()
  centerOfMassCoord = 0.d0

  CALL faceEntityIt%next(faceEntityKey, faceEntityValue)

  element = faceEntityKey%ints(1)

  do localNodeId = 1, 3
    globalNodeId = globalNodeIdTypeFace(localNodeId,element)

    centerOfMassCoord = centerOfMassCoord + nodeCoord(axis,globalNodeId)
  enddo

  centerOfMassCoord = centerOfMassCoord / 3.0d0

  if (ABS(centerOfMassCoord - boundaryCoord) < tol) then
    faceKey%ints(1) = element

    CALL faceHash%set(faceKey, faceEntityValue)
    cycle
  endif
enddo

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshFaceEntities
