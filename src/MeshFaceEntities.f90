!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshFaceEntities(axis, face_entity_hash, face1_hash, face2_hash)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: periodicFaceId
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
character(len=1), intent(in) :: axis

type(fhash_type__ints_double), intent(inout) :: face_entity_hash
type(fhash_type_iterator__ints_double)       :: face_entity_it
type(ints_type)                              :: face_entity_key
integer                                      :: face_entity_value

type(fhash_type__ints_double), intent(out) :: face1_hash, face2_hash
type(ints_type)                            :: face1_key, face2_key
integer                                    :: face1_size, face2_size

integer :: kk, periodic_face1, periodic_face2
!----------------------------------------------------------------------------------------------------------------------------------!
face1_size = 0
face2_size = 0

periodic_face1 = 0
periodic_face2 = 0

if (axis=='x') then
  periodic_face1 = periodicFaceId(1)
  periodic_face2 = periodicFaceId(2)
elseif (axis=='y') then
  periodic_face1 = periodicFaceId(3)
  periodic_face2 = periodicFaceId(4)
elseif (axis=='z') then
  periodic_face1 = periodicFaceId(5)
  periodic_face2 = periodicFaceId(6)
endif

! Calculate the size of face entity hashes
call face_entity_it%begin(face_entity_hash)
do kk = 1, face_entity_hash%key_count()
  call face_entity_it%next(face_entity_key, face_entity_value)

  if (face_entity_value==periodic_face1) face1_size = face1_size + 1
  if (face_entity_value==periodic_face2) face2_size = face2_size + 1
enddo

allocate(face1_key%ints(1))
allocate(face2_key%ints(1))

call face1_hash%reserve(face1_size)
call face2_hash%reserve(face2_size)

! Fill the face entity hashes
call face_entity_it%begin(face_entity_hash)
do kk = 1, face_entity_hash%key_count()
  call face_entity_it%next(face_entity_key, face_entity_value)
  if (face_entity_value==periodic_face1) then
    face1_key%ints(1) = face_entity_key%ints(1)
    call face1_hash%set(face1_key, face_entity_value)
    cycle
  endif
  if (face_entity_value==periodic_face2) then
    face2_key%ints(1) = face_entity_key%ints(1)
    call face2_hash%set(face2_key, face_entity_value)
    cycle
  endif
enddo

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshFaceEntities
