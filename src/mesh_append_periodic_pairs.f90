!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine mesh_append_periodic_pairs(elemcon, starting_pair, ending_pair, node_pairing_hash)
!----------------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use kcw_mod, only: F_m
use geometry_mod, only: node_pair_id
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)               :: elemcon_key
integer                       :: elemcon_value

integer, intent(in) :: starting_pair, ending_pair

type(fhash_type__ints_double), intent(inout) :: node_pairing_hash
type(fhash_type_iterator__ints_double)     :: node_pairing_it
type(ints_type)                            :: node_pairing_key
integer                                    :: node_pairing_value

integer :: kk, pp
integer :: source, dest

logical :: success
!----------------------------------------------------------------------------------------------------------------------------------!
allocate(elemcon_key%ints(2))

call node_pairing_it%begin(node_pairing_hash)

do kk = starting_pair, ending_pair
    call node_pairing_it%next(node_pairing_key, node_pairing_value)

    source = node_pairing_key%ints(1)
    dest   = node_pairing_value

    ! Append pair
    F_m%row(kk) = source
    F_m%col(kk) = dest

    elemcon_key%ints(1) = source
    elemcon_key%ints(2) = dest

    call elemcon%get(elemcon_key, elemcon_value, success)

    if (success) then
        node_pair_id(kk) = elemcon_value
    else
        call elemcon%set(elemcon_key, kk)
        node_pair_id(kk) = kk
    endif

    ! Append inverse pair
    pp = kk + node_pairing_hash%key_count()

    F_m%row(pp) = dest
    F_m%col(pp) = source

    elemcon_key%ints(1) = dest
    elemcon_key%ints(2) = source

    call elemcon%get(elemcon_key, elemcon_value, success)

    if (success) then
        node_pair_id(pp) = elemcon_value
    else
        call elemcon%set(elemcon_key, pp)
        node_pair_id(pp) = pp
    endif
enddo

deallocate(elemcon_key%ints)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine mesh_append_periodic_pairs
