!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshDirichletFaces(numElementsTypeFace, nenTypeFace, globalNodeIdTypeFace, faceEntityHash)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: numDirichletFaces, numNanopFaces, dirichletFaceId, nanopFaceId
use geometry_mod,    only: numDimensions, nodeBelongsToDirichletFace, nodeBelongsToFaceId, numNodes, &
                           isDirichletFace, boxLow, boxHigh, nodeCoord
use iofiles_mod,     only: IO_dirichletFaces
use constants_mod,   only: tol
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
integer, intent(in)                                              :: numElementsTypeFace, nenTypeFace
integer, dimension(nenTypeFace, numElementsTypeFace), intent(in) :: globalNodeIdTypeFace
integer                                                          :: face, element, axis, counter, localNodeId, globalNodeId

type(fhash_type__ints_double), intent(inout) :: faceEntityHash
type(ints_type)                              :: faceEntityKey
integer                                      :: faceEntityValue
!----------------------------------------------------------------------------------------------------------------------------------!
allocate(nodeBelongsToDirichletFace(numNodes))
allocate(nodeBelongsToFaceId(numNodes))

nodeBelongsToDirichletFace = .False.
nodeBelongsToFaceId        = -1

allocate(faceEntityKey%ints(1))

#ifdef DEBUG_OUTPUTS
open(unit=123, file = IO_dirichletFaces)
#endif

isDirichletFace = .False.

do element = 1, numElementsTypeFace
  faceEntityKey%ints(1) = element
  call faceEntityHash%get(faceEntityKey, faceEntityValue)
  do face = 1, numDirichletFaces
    if (faceEntityValue == dirichletFaceId(face)) then
      do localNodeId = 1, nenTypeFace
        globalNodeId = globalNodeIdTypeFace(localNodeId,element)

        nodeBelongsToDirichletFace(globalNodeId) = .True.
        nodeBelongsToFaceId(globalNodeId)        = faceEntityValue

        ! Find if a node is located at a corner
        counter = 0
        do axis = 1, numDimensions
          if (DABS(nodeCoord(axis,globalNodeId) - boxLow(axis))  < tol) counter = counter + 1
          if (DABS(nodeCoord(axis,globalNodeId) - boxHigh(axis)) < tol) counter = counter + 1
        enddo

        ! If a node is located at a corner skip the loop
        if (counter > 1) cycle

        do axis = 1, numDimensions
          if (DABS(nodeCoord(axis,globalNodeId) - boxLow(axis))  < tol) isDirichletFace(axis,1) = .True.
          if (DABS(nodeCoord(axis,globalNodeId) - boxHigh(axis)) < tol) isDirichletFace(axis,2) = .True.
        enddo
#ifdef DEBUG_OUTPUTS
        write(123,'(3(E16.8),6(L3))') nodeCoord(1,globalNodeId), nodeCoord(2,globalNodeId), nodeCoord(3,globalNodeId), &
                                      (isDirichletFace(axis,1), axis = 1, 3), (isDirichletFace(axis,2), axis = 1, 3)
#endif
      enddo
    endif
  enddo

  ! Nanops section
  do face = 1, numNanopFaces
    if (faceEntityValue == nanopFaceId(face)) then
      do localNodeId = 1, nenTypeFace
        globalNodeId = globalNodeIdTypeFace(localNodeId,element)

        nodeBelongsToDirichletFace(globalNodeId) = .True.
        nodeBelongsToFaceId(globalNodeId)        = faceEntityValue
      enddo
    endif
  enddo
enddo

#ifdef DEBUG_OUTPUTS
close(123)
#endif

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshDirichletFaces
