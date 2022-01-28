subroutine mesh_append_dest_neighbors(elemcon, forward_steps, num_dest_neighbors, node_pairing_hash)

use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use kcw_mod,      only: F_m
use geometry_mod, only: node_pair_id, num_of_bulk_pairs

implicit none

type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemcon_key
integer                                      :: elemcon_value

integer, intent(in) :: forward_steps, num_dest_neighbors

type(fhash_type__ints_double), intent(inout) :: node_pairing_hash
type(fhash_type_iterator__ints_double)       :: node_pairing_it
type(ints_type)                              :: node_pairing_key
integer                                      :: node_pairing_value

integer :: kk, pp, node_pair, aux
integer :: source, dest

logical :: success

aux = 0

allocate(elemcon_key%ints(2))

call node_pairing_it%begin(node_pairing_hash)

do kk = 1, node_pairing_hash%key_count()
    call node_pairing_it%next(node_pairing_key, node_pairing_value)

    source = node_pairing_key%ints(1)
    dest   = node_pairing_value

    do node_pair = 1, num_of_bulk_pairs + forward_steps
        if ((F_m%col(node_pair)==dest).and.(F_m%row(node_pair).ne.dest)) then

            if (F_m%is_zero(node_pair)) cycle

            aux = aux + 1

            pp = num_of_bulk_pairs + forward_steps + 2*node_pairing_hash%key_count() + aux

            ! Append pair
            F_m%row(pp) = source
            F_m%col(pp) = F_m%row(node_pair)

            elemcon_key%ints(1) = source
            elemcon_key%ints(2) = F_m%row(node_pair)

            call elemcon%get(elemcon_key, elemcon_value, success)

            if (success) then
                node_pair_id(pp) = elemcon_value
            else
                call elemcon%set(elemcon_key, pp)
                node_pair_id(pp) = pp
            endif

            ! Append inverse pair
            F_m%row(pp + num_dest_neighbors) = F_m%row(node_pair)
            F_m%col(pp + num_dest_neighbors) = source

            elemcon_key%ints(1) = F_m%row(node_pair)
            elemcon_key%ints(2) = source

            call elemcon%get(elemcon_key, elemcon_value, success)

            if (success) then
                node_pair_id(pp + num_dest_neighbors) = elemcon_value
            else
                call elemcon%set(elemcon_key, pp + num_dest_neighbors)
                node_pair_id(pp + num_dest_neighbors) = pp + num_dest_neighbors
            endif
        endif
    enddo
enddo

deallocate(elemcon_key%ints)

return
end subroutine mesh_append_dest_neighbors
