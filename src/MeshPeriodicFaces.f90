!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshPeriodicFaces(axis, axisFaceOneHash, axisFaceTwoHash)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use iofiles_mod, only: IO_xFaceOneElements, IO_xFaceTwoElements, &
                       IO_yFaceOneElements, IO_yFaceTwoElements, &
                       IO_zFaceOneElements, IO_zFaceTwoElements
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
character(len=1), intent(in) :: axis

type(fhash_type__ints_double), intent(inout) :: axisFaceOneHash, axisFaceTwoHash
type(fhash_type_iterator__ints_double)       :: axisFaceOneIt, axisFaceTwoIt
type(ints_type)                              :: axisFaceOneKey, axisFaceTwoKey
integer                                      :: axisFaceOneValue, axisFaceTwoValue

integer :: keyIndex
!----------------------------------------------------------------------------------------------------------------------------------!
if (axis=='x') then
  open(unit=111, file=IO_xFaceOneElements)
  open(unit=222, file=IO_xFaceTwoElements)
elseif (axis=='y') then
  open(unit=111, file=IO_yFaceOneElements)
  open(unit=222, file=IO_yFaceTwoElements)
elseif (axis=='z') then
  open(unit=111, file=IO_zFaceOneElements)
  open(unit=222, file=IO_zFaceTwoElements)
endif

call axisFaceOneIt%begin(axisFaceOneHash)
call axisFaceTwoIt%begin(axisFaceTwoHash)

do keyIndex = 1, axisFaceOneHash%key_count()
  call axisFaceOneIt%next(axisFaceOneKey, axisFaceOneValue)
  write(111,*) axisFaceOneValue, axisFaceOneKey%ints(1)
enddo

do keyIndex = 1, axisFaceTwoHash%key_count()
  call axisFaceTwoIt%next(axisFaceTwoKey, axisFaceTwoValue)
  write(222,*) axisFaceTwoValue, axisFaceTwoKey%ints(1)
enddo

close(111)
close(222)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicFaces
