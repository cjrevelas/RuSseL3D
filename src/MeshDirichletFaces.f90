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
integer                                                          :: ii, jj, kk, mm, pp, idummy

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

do jj = 1, numElementsTypeFace
  faceEntityKey%ints(1) = jj
  call faceEntityHash%get(faceEntityKey, faceEntityValue)
  do ii = 1, numDirichletFaces
    if (faceEntityValue == dirichletFaceId(ii)) then
      do pp = 1, nenTypeFace
        idummy = globalNodeIdTypeFace(pp,jj)

        nodeBelongsToDirichletFace(idummy) = .True.
        nodeBelongsToFaceId(idummy)        = faceEntityValue

        ! Find if a node is located at a corner
        kk = 0
        do mm = 1, numDimensions
          if (DABS(nodeCoord(mm,idummy) - boxLow(mm))  < tol) kk = kk + 1
          if (DABS(nodeCoord(mm,idummy) - boxHigh(mm)) < tol) kk = kk + 1
        enddo

        ! If a node is located at a corner skip the loop
        if (kk > 1) cycle

        do mm = 1, numDimensions
          if (DABS(nodeCoord(mm,idummy) - boxLow(mm))  < tol) isDirichletFace(mm,1) = .True.
          if (DABS(nodeCoord(mm,idummy) - boxHigh(mm)) < tol) isDirichletFace(mm,2) = .True.
        enddo
#ifdef DEBUG_OUTPUTS
        write(123,'(3(E16.8),6(L3))') nodeCoord(1,idummy), nodeCoord(2,idummy), nodeCoord(3,idummy), &
                                      (isDirichletFace(kk,1), kk = 1,3), (isDirichletFace(kk,2), kk = 1,3)
#endif
      enddo
    endif
  enddo

  ! Nanops section
  do ii = 1, numNanopFaces
    if (faceEntityValue == nanopFaceId(ii)) then
      do kk = 1, nenTypeFace
        idummy = globalNodeIdTypeFace(kk,jj)

        nodeBelongsToDirichletFace(idummy) = .True.
        nodeBelongsToFaceId(idummy)        = faceEntityValue
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
