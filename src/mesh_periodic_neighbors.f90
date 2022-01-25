subroutine mesh_periodic_neighbors(node_pairing_hash, num_dest_neighbors)

use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module

use geometry_mod, only: num_of_elems_of_node

implicit none

type(fhash_type__ints_double), intent(inout) :: node_pairing_hash
type(fhash_type_iterator__ints_double)       :: node_pairing_it
type(ints_type)                              :: node_pairing_key
integer                                      :: node_pairing_value

integer, intent(out) :: num_dest_neighbors

integer :: kk, source_xx, dest_xx

call node_pairing_it%begin(node_pairing_hash)

do kk = 1, node_pairing_hash%key_count()
    call node_pairing_it%next(node_pairing_key, node_pairing_value)

    source_xx = node_pairing_key%ints(1)
    dest_xx   = node_pairing_value

    num_dest_neighbors = num_dest_neighbors + 3 * num_of_elems_of_node(dest_xx)
enddo

return
end subroutine mesh_periodic_neighbors
