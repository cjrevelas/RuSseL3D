!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine mesh_dirichlet_faces(numel_type_face, nen_type_face, global_node_id_type_face, face_entity_hash)
!----------------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: numDirichletFaces, numNanoparticleFaces, dirichletFaceId, nanoparticleFaceId
use geometry_mod,    only: ndm, nodeBelongsToDirichletFace, numnp, isDirichletFace, boxLow, boxHigh, xc
use iofiles_mod,     only: dir_faces
use constants_mod,   only: tol
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numel_type_face, nen_type_face
integer, dimension(nen_type_face, numel_type_face), intent(in) :: global_node_id_type_face
integer :: ii, jj, kk, mm, pp, idummy

type(fhash_type__ints_double), intent(inout) :: face_entity_hash
type(ints_type) :: face_entity_key
integer :: face_entity_value
!----------------------------------------------------------------------------------------------------------------------------------!
allocate(nodeBelongsToDirichletFace(numnp))
nodeBelongsToDirichletFace = .False.

allocate(face_entity_key%ints(1))

#ifdef DEBUG_OUTPUTS
open(unit=123, file = dir_faces)
#endif

isDirichletFace = .False.

do jj = 1, numel_type_face
  face_entity_key%ints(1) = jj
  call face_entity_hash%get(face_entity_key, face_entity_value)
  do ii = 1, numDirichletFaces
    if (face_entity_value == dirichletFaceId(ii)) then
      do pp = 1, nen_type_face
        idummy = global_node_id_type_face(pp,jj)

        nodeBelongsToDirichletFace(idummy) = .True.

        ! Find if a node is located at a corner
        kk = 0
        do mm = 1, ndm
          if (DABS(xc(mm, idummy) - boxLow(mm)) < tol) kk = kk + 1
          if (DABS(xc(mm, idummy) - boxHigh(mm)) < tol) kk = kk + 1
        enddo

        ! If a node is located at a corner skip the loop
        if (kk > 1) cycle

        do mm = 1, ndm
          if (DABS(xc(mm, idummy) - boxLow(mm)) < tol) isDirichletFace(mm,1) = .True.
          if (DABS(xc(mm, idummy) - boxHigh(mm)) < tol) isDirichletFace(mm,2) = .True.
        enddo
#ifdef DEBUG_OUTPUTS
        write(123,'(3(E16.8),6(L3))') xc(1, idummy), xc(2, idummy), xc(3, idummy), (isDirichletFace(kk,1), kk = 1,3), (isDirichletFace(kk,2), kk = 1,3)
#endif
      enddo
    endif
  enddo

  ! Nanoparticles section
  do ii = 1, numNanoparticleFaces
    if (face_entity_value==nanoparticleFaceId(ii)) then
      do kk = 1, nen_type_face
        idummy = global_node_id_type_face(kk,jj)

        nodeBelongsToDirichletFace(idummy) = .True.
      enddo
    endif
  enddo
enddo

#ifdef DEBUG_OUTPUTS
close(123)
#endif

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine mesh_dirichlet_faces
