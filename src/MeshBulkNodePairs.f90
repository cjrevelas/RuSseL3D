!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshBulkNodePairs(nodePairingXXhash, nodePairingYYhash, nodePairingZZhash, elemcon)
!----------------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: periodicAxisId
use geometry_mod,    only: numNodesLocalTypeDomain, numElementsTypeDomain, numTotalNodePairs, nodePairId, globalNodeIdTypeDomain
use kcw_mod,         only: F_m
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
logical :: success

type(fhash_type__ints_double), intent(inout) :: nodePairingXXhash, nodePairingYYhash, nodePairingZZhash, elemcon
type(fhash_type_iterator__ints_double)       :: nodePairingXXit, nodePairingYYit, nodePairingZZit
type(ints_type)                              :: nodePairingXXkey, nodePairingYYkey, nodePairingZZkey, elemcon_key
integer                                      :: nodePairingXXvalue, nodePairingYYvalue, nodePairingZZvalue, elemcon_num_of_keys, elemcon_value

integer :: ii, jj, mm, kk
integer :: node_pair
integer :: source_xx, dest_xx
integer :: source_yy, dest_yy
integer :: source_zz, dest_zz
!----------------------------------------------------------------------------------------------------------------------------------!
allocate(nodePairId(numTotalNodePairs))
nodePairId = 0

! Assembly the nodePairId (hash) matrix
allocate(elemcon_key%ints(2)) ! Each key is defined by a pair (2) of nodes

! Total number of required keys
elemcon_num_of_keys = 2 * numNodesLocalTypeDomain * numElementsTypeDomain
call elemcon%reserve(elemcon_num_of_keys)

node_pair = 0
do mm = 1, numElementsTypeDomain
  do jj = 1, numNodesLocalTypeDomain
    do ii = 1, numNodesLocalTypeDomain
      node_pair = node_pair + 1

      ! Define the pair of nodes to be examined and assigned a elemcon_value
      elemcon_key%ints(1) = globalNodeIdTypeDomain(jj,mm)
      elemcon_key%ints(2) = globalNodeIdTypeDomain(ii,mm)

      if (periodicAxisId(1)) then
        call nodePairingXXit%begin(nodePairingXXhash)
        do kk = 1, nodePairingXXhash%key_count()
          call nodePairingXXit%next(nodePairingXXkey, nodePairingXXvalue)
          source_xx = nodePairingXXkey%ints(1)
          dest_xx   = nodePairingXXvalue

          if ((elemcon_key%ints(1) == dest_xx).AND.(elemcon_key%ints(2) == dest_xx))   exit
          if ((elemcon_key%ints(1) == dest_xx).AND.(elemcon_key%ints(2) == source_xx)) exit

          if (elemcon_key%ints(1) == dest_xx) elemcon_key%ints(1) = source_xx
          if (elemcon_key%ints(2) == dest_xx) elemcon_key%ints(2) = source_xx
        enddo
      endif

      if (periodicAxisId(2)) then
        call nodePairingYYit%begin(nodePairingYYhash)
        do kk = 1, nodePairingYYhash%key_count()
          call nodePairingYYit%next(nodePairingYYkey, nodePairingYYvalue)
          source_yy = nodePairingYYkey%ints(1)
          dest_yy   = nodePairingYYvalue

          if ((elemcon_key%ints(1) == dest_yy).AND.(elemcon_key%ints(2) == dest_yy))   exit
          if ((elemcon_key%ints(1) == dest_yy).AND.(elemcon_key%ints(2) == source_yy)) exit

          if (elemcon_key%ints(1) == dest_yy) elemcon_key%ints(1) = source_yy
          if (elemcon_key%ints(2) == dest_yy) elemcon_key%ints(2) = source_yy
        enddo
      endif

      if (periodicAxisId(3)) then
        call nodePairingZZit%begin(nodePairingZZhash)
        do kk = 1, nodePairingZZhash%key_count()
          call nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)
          source_zz = nodePairingZZkey%ints(1)
          dest_zz   = nodePairingZZvalue

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
        nodePairId(node_pair) = elemcon_value  ! This pair has already been met, thus assigned a elemcon_value
      else
        call elemcon%set(elemcon_key, node_pair) ! Store the new elemcon_value for next iteration's check
        nodePairId(node_pair) = node_pair        ! This pair is met for the first time
      endif
    enddo
  enddo
enddo

deallocate(elemcon_key%ints)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshBulkNodePairs
