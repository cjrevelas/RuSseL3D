!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine mesh_bulk_node_pairs(elemcon)
!----------------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use geometry_mod, only: nel, numel, total_num_of_node_pairs, node_pair_id, global_node_id_type_domain
use kcw_mod,      only: F_m
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
logical :: success

type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)               :: elemcon_key
integer                       :: elemcon_num_of_keys, elemcon_value

integer :: ii, jj, mm
integer :: node_pair
!----------------------------------------------------------------------------------------------------------------------------------!
allocate(node_pair_id(total_num_of_node_pairs))
node_pair_id = 0

! Assembly the node_pair_id (hash) matrix
allocate(elemcon_key%ints(2)) ! Each key is defined by a pair (2) of nodes

! Total number of required keys
elemcon_num_of_keys = 2 * nel * numel
call elemcon%reserve(elemcon_num_of_keys)

node_pair = 0
do mm = 1, numel
    do jj = 1, nel
        do ii = 1, nel
            node_pair = node_pair + 1

            F_m%row(node_pair) = global_node_id_type_domain(jj,mm)
            F_m%col(node_pair) = global_node_id_type_domain(ii,mm)

            ! Define the pair of nodes to be examined and assigned a elemcon_value
            elemcon_key%ints(1) = global_node_id_type_domain(jj,mm)
            elemcon_key%ints(2) = global_node_id_type_domain(ii,mm)

            ! Assign elemcon_value to the pair
            call elemcon%get(elemcon_key, elemcon_value, success)

            if (success) then
               node_pair_id(node_pair) = elemcon_value  ! This pair has already been met, thus assigned a elemcon_value
            else
               call elemcon%set(elemcon_key, node_pair) ! Store the new elemcon_value for next iteration's check
               node_pair_id(node_pair) = node_pair      ! This pair is met for the first time
            endif
        enddo
    enddo
enddo

deallocate(elemcon_key%ints)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine mesh_bulk_node_pairs
