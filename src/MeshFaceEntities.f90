!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshFaceEntities(axis, faceEntityHash, faceOneHash, faceTwoHash)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: periodicFaceId
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
character(len=1), intent(in) :: axis

type(fhash_type__ints_double), intent(inout) :: faceEntityHash
type(fhash_type_iterator__ints_double)       :: faceEntityIt
type(ints_type)                              :: faceEntityKey
integer                                      :: faceEntityValue

type(fhash_type__ints_double), intent(out) :: faceOneHash, faceTwoHash
type(ints_type)                            :: faceOneKey, faceTwoKey
integer                                    :: faceOneSize, faceTwoSize

integer :: kk, periodicFaceOne, periodicFaceTwo
!----------------------------------------------------------------------------------------------------------------------------------!
faceOneSize = 0
faceTwoSize = 0

periodicFaceOne = 0
periodicFaceTwo = 0

if (axis=='x') then
  periodicFaceOne = periodicFaceId(1)
  periodicFaceTwo = periodicFaceId(2)
elseif (axis=='y') then
  periodicFaceOne = periodicFaceId(3)
  periodicFaceTwo = periodicFaceId(4)
elseif (axis=='z') then
  periodicFaceOne = periodicFaceId(5)
  periodicFaceTwo = periodicFaceId(6)
endif

! Calculate the size of face entity hashes
call faceEntityIt%begin(faceEntityHash)
do kk = 1, faceEntityHash%key_count()
  call faceEntityIt%next(faceEntityKey, faceEntityValue)

  if (faceEntityValue == periodicFaceOne) faceOneSize = faceOneSize + 1
  if (faceEntityValue == periodicFaceTwo) faceTwoSize = faceTwoSize + 1
enddo

allocate(faceOneKey%ints(1))
allocate(faceTwoKey%ints(1))

call faceOneHash%reserve(faceOneSize)
call faceTwoHash%reserve(faceTwoSize)

! Fill the face entity hashes
call faceEntityIt%begin(faceEntityHash)
do kk = 1, faceEntityHash%key_count()
  call faceEntityIt%next(faceEntityKey, faceEntityValue)
  if (faceEntityValue == periodicFaceOne) then
    faceOneKey%ints(1) = faceEntityKey%ints(1)
    call faceOneHash%set(faceOneKey, faceEntityValue)
    cycle
  endif
  if (faceEntityValue == periodicFaceTwo) then
    faceTwoKey%ints(1) = faceEntityKey%ints(1)
    call faceTwoHash%set(faceTwoKey, faceEntityValue)
    cycle
  endif
enddo

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshFaceEntities
