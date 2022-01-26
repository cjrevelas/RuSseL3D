subroutine mesh_build_node_pairing(global_node_id_type_face, axis, face1_hash, face2_hash, node_pairing_hash)

use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module

use geometry_mod, only: xc, nen_type_face, numel_type_face
use iofiles_mod,  only: node_pairing_xx, node_pairing_yy, node_pairing_zz

implicit none

character(len=1), intent(in) :: axis

integer, dimension(nen_type_face, numel_type_face), intent(in) :: global_node_id_type_face

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
integer :: dimension1, dimension2

if (axis=='x') then
    dimension1 = 2 ! y
    dimension2 = 3 ! z

    open(unit=1111, file=node_pairing_xx)
elseif (axis=='y') then
    dimension1 = 1 ! x
    dimension2 = 3 ! z

    open(unit=1111, file=node_pairing_yy)
elseif (axis=='z') then
    dimension1 = 1 ! x
    dimension2 = 2 !y

    open(unit=1111, file=node_pairing_zz)
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
            node1 = global_node_id_type_face(ii,elem1)
            do jj = 1, 3
                node2 = global_node_id_type_face(jj, elem2)
                if ((ABS(xc(dimension1,node1)-xc(dimension1,node2))<1e-12).AND.(ABS(xc(dimension2,node1)-xc(dimension2,node2))<1e-12)) then
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
end subroutine mesh_build_node_pairing
