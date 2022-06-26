!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine FemApplyPeriodicity(node_pairing_hash, elemcon)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod, only: F_m
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: node_pairing_hash
type(fhash_type_iterator__ints_double)       :: node_pairing_it
type(ints_type)                              :: node_pairing_key
integer                                      :: node_pairing_value

type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemcon_key

integer :: source, dest
integer :: kk, mm1, mm2, mm3
!----------------------------------------------------------------------------------------------------------------------------------!
call node_pairing_it%begin(node_pairing_hash)

allocate(elemcon_key%ints(2))

do kk = 1, node_pairing_hash%key_count()
  call node_pairing_it%next(node_pairing_key, node_pairing_value)

  source = node_pairing_key%ints(1)
  dest   = node_pairing_value

  elemcon_key%ints(1) = source
  elemcon_key%ints(2) = source
  call elemcon%get(elemcon_key, mm1)

  elemcon_key%ints(1) = dest
  elemcon_key%ints(2) = dest
  call elemcon%get(elemcon_key, mm2)

  F_m%g(mm1)  = F_m%g(mm1)  + F_m%g(mm2)
  F_m%rh(mm1) = F_m%rh(mm1) + F_m%rh(mm2)
  F_m%g(mm2)  = 1.0d0
  F_m%rh(mm2) = 0.0d0

  elemcon_key%ints(1) = dest
  elemcon_key%ints(2) = source
  call elemcon%get(elemcon_key, mm3)

  F_m%g(mm3) = -1.0d0
enddo

deallocate(elemcon_key%ints)
return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine FemApplyPeriodicity
