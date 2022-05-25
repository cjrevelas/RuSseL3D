!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine mesh_bulk_node_pairs(node_pairing_xx_hash, node_pairing_yy_hash, node_pairing_zz_hash, elemcon)
!----------------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: periodicAxisId
use geometry_mod, only: nel, numel, numTotalNodePairs, node_pair_id, global_node_id_type_domain
use kcw_mod,      only: F_m
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
logical :: success

type(fhash_type__ints_double), intent(inout) :: node_pairing_xx_hash, node_pairing_yy_hash, node_pairing_zz_hash, elemcon
type(fhash_type_iterator__ints_double)       :: node_pairing_xx_it, node_pairing_yy_it, node_pairing_zz_it
type(ints_type)                              :: node_pairing_xx_key, node_pairing_yy_key, node_pairing_zz_key, elemcon_key
integer                                      :: node_pairing_xx_value, node_pairing_yy_value, node_pairing_zz_value, &
&                                               elemcon_num_of_keys, elemcon_value

integer :: ii, jj, mm, kk
integer :: node_pair
integer :: source_xx, dest_xx
integer :: source_yy, dest_yy
integer :: source_zz, dest_zz
!----------------------------------------------------------------------------------------------------------------------------------!
allocate(node_pair_id(numTotalNodePairs))
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

      ! Define the pair of nodes to be examined and assigned a elemcon_value
      elemcon_key%ints(1) = global_node_id_type_domain(jj,mm)
      elemcon_key%ints(2) = global_node_id_type_domain(ii,mm)

      if (periodicAxisId(1)) then
        call node_pairing_xx_it%begin(node_pairing_xx_hash)
        do kk = 1, node_pairing_xx_hash%key_count()
          call node_pairing_xx_it%next(node_pairing_xx_key, node_pairing_xx_value)
          source_xx = node_pairing_xx_key%ints(1)
          dest_xx   = node_pairing_xx_value

          if ((elemcon_key%ints(1) == dest_xx).AND.(elemcon_key%ints(2) == dest_xx))   exit
          if ((elemcon_key%ints(1) == dest_xx).AND.(elemcon_key%ints(2) == source_xx)) exit

          if (elemcon_key%ints(1) == dest_xx) elemcon_key%ints(1) = source_xx
          if (elemcon_key%ints(2) == dest_xx) elemcon_key%ints(2) = source_xx
        enddo
      endif

      if (periodicAxisId(2)) then
        call node_pairing_yy_it%begin(node_pairing_yy_hash)
        do kk = 1, node_pairing_yy_hash%key_count()
          call node_pairing_yy_it%next(node_pairing_yy_key, node_pairing_yy_value)
          source_yy = node_pairing_yy_key%ints(1)
          dest_yy   = node_pairing_yy_value

          if ((elemcon_key%ints(1) == dest_yy).AND.(elemcon_key%ints(2) == dest_yy))   exit
          if ((elemcon_key%ints(1) == dest_yy).AND.(elemcon_key%ints(2) == source_yy)) exit

          if (elemcon_key%ints(1) == dest_yy) elemcon_key%ints(1) = source_yy
          if (elemcon_key%ints(2) == dest_yy) elemcon_key%ints(2) = source_yy
        enddo
      endif

      if (periodicAxisId(3)) then
        call node_pairing_zz_it%begin(node_pairing_zz_hash)
        do kk = 1, node_pairing_zz_hash%key_count()
          call node_pairing_zz_it%next(node_pairing_zz_key, node_pairing_zz_value)
          source_zz = node_pairing_zz_key%ints(1)
          dest_zz   = node_pairing_zz_value

          if ((elemcon_key%ints(1) == dest_zz).AND.(elemcon_key%ints(2) == dest_zz))   exit
          if ((elemcon_key%ints(1) == dest_zz).AND.(elemcon_key%ints(2) == source_zz)) exit

          if (elemcon_key%ints(1) == dest_zz) elemcon_key%ints(1) = source_zz
          if (elemcon_key%ints(2) == dest_zz) elemcon_key%ints(2) = source_zz
        enddo
      endif

      F_m%row(node_pair) = elemcon_key%ints(1)
      F_m%col(node_pair) = elemcon_key%ints(2)

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
