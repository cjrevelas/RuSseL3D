!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshPeriodicFaces(axis, axis_face1_hash, axis_face2_hash)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use iofiles_mod, only: IO_xFace1Elements, IO_xFace2Elements, &
                       IO_yFace1Elements, IO_yFace2Elements, &
                       IO_zFace1Elements, IO_zFace2Elements
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
character(len=1), intent(in) :: axis

type(fhash_type__ints_double), intent(inout) :: axis_face1_hash, axis_face2_hash
type(fhash_type_iterator__ints_double)       :: axis_face1_it, axis_face2_it
type(ints_type)                              :: axis_face1_key, axis_face2_key
integer                                      :: axis_face1_value, axis_face2_value

integer :: kk
!----------------------------------------------------------------------------------------------------------------------------------!
if (axis=='x') then
  open(unit=111, file=IO_xFace1Elements)
  open(unit=222, file=IO_xFace2Elements)
elseif (axis=='y') then
  open(unit=111, file=IO_yFace1Elements)
  open(unit=222, file=IO_yFace2Elements)
elseif (axis=='z') then
  open(unit=111, file=IO_zFace1Elements)
  open(unit=222, file=IO_zFace2Elements)
endif

call axis_face1_it%begin(axis_face1_hash)
call axis_face2_it%begin(axis_face2_hash)

do kk = 1, axis_face1_hash%key_count()
  call axis_face1_it%next(axis_face1_key, axis_face1_value)
  write(111,*) axis_face1_value, axis_face1_key%ints(1)
enddo

do kk = 1, axis_face2_hash%key_count()
  call axis_face2_it%next(axis_face2_key, axis_face2_value)
  write(222,*) axis_face2_value, axis_face2_key%ints(1)
enddo

close(111)
close(222)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicFaces
