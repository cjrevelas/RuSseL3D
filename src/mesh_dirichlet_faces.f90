subroutine mesh_dirichlet_faces(numel_type_face, nen_type_face, global_node_id_type_face, face_entity_hash)

use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: num_of_dirichlet_faces, num_of_nanoparticle_faces, ids_dirichlet_faces, ids_nanopart_faces
use geometry_mod,    only: ndm, node_belongs_to_dirichlet_face, numnp, is_dirichlet_face, box_lo, box_hi, xc
use iofiles_mod,     only: dir_faces
use constants_mod,   only: tol

implicit none

integer, intent(in) :: numel_type_face, nen_type_face
integer, dimension(nen_type_face, numel_type_face), intent(in) :: global_node_id_type_face
integer :: ii, jj, kk, mm, pp, idummy

type(fhash_type__ints_double), intent(inout) :: face_entity_hash
type(ints_type) :: face_entity_key
integer :: face_entity_value

allocate(node_belongs_to_dirichlet_face(numnp))
node_belongs_to_dirichlet_face = .false.

allocate(face_entity_key%ints(1))

#ifdef DEBUG_OUTPUTS
open(unit=123, file = dir_faces)
#endif

is_dirichlet_face = .false.

do jj = 1, numel_type_face
    face_entity_key%ints(1) = jj
    call face_entity_hash%get(face_entity_key, face_entity_value)
    do ii = 1, num_of_dirichlet_faces
        if (face_entity_value==ids_dirichlet_faces(ii)) then
            do pp = 1, nen_type_face
                idummy = global_node_id_type_face(pp,jj)

                node_belongs_to_dirichlet_face(idummy) = .True.

                ! Find if a node is located at a corner
                kk = 0
                do mm = 1, ndm
                    if (DABS(xc(mm, idummy) - box_lo(mm)) < tol) then
                        kk = kk + 1
                    endif
                    if (DABS(xc(mm, idummy) - box_hi(mm)) < tol) then
                        kk = kk + 1
                    endif
                enddo

                ! If a node is located at a corner skip the loop
                if (kk > 1) cycle

                do mm = 1, ndm
                    if (DABS(xc(mm, idummy) - box_lo(mm)) < tol) then
                        is_dirichlet_face(mm,1) = .true.
                    endif
                    if (DABS(xc(mm, idummy) - box_hi(mm)) < tol) then
                        is_dirichlet_face(mm,2) = .true.
                    endif
                enddo
#ifdef DEBUG_OUTPUTS
                write(123,'(3(E16.8),6(L3))') xc(1, idummy), xc(2, idummy), xc(3, idummy), &
    &                                 (is_dirichlet_face(kk,1), kk = 1,3), (is_dirichlet_face(kk,2), kk = 1,3)
#endif
            enddo
        endif
    enddo

    ! Nanoparticles section
    do ii = 1, num_of_nanoparticle_faces
        if (face_entity_value==ids_nanopart_faces(ii)) then
            do kk = 1, nen_type_face
                idummy = global_node_id_type_face(kk,jj)

                node_belongs_to_dirichlet_face(idummy) = .True.
            enddo
        endif
    enddo
enddo

#ifdef DEBUG_OUTPUTS
close(123)
#endif

return
end subroutine mesh_dirichlet_faces
