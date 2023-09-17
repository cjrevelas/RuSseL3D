!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine FemApplyPeriodicity(nodePairingHash, elemcon)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod, only: F_m
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: nodePairingHash
type(fhash_type_iterator__ints_double)       :: nodePairingIt
type(ints_type)                              :: nodePairingKey
integer                                      :: nodePairingValue

type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemconKey

integer :: src, dst
integer :: keyIndex, pairOne, pairTwo, pairThree
!----------------------------------------------------------------------------------------------------------------------------------!
CALL nodePairingIt%begin(nodePairingHash)

allocate(elemconKey%ints(2))

do keyIndex = 1, nodePairingHash%key_count()
  CALL nodePairingIt%next(nodePairingKey, nodePairingValue)

  src = nodePairingKey%ints(1)
  dst = nodePairingValue

  elemconKey%ints(1) = src
  elemconKey%ints(2) = src
  CALL elemcon%get(elemconKey, pairOne)

  elemconKey%ints(1) = dst
  elemconKey%ints(2) = dst
  CALL elemcon%get(elemconKey, pairTwo)

  F_m%g(pairOne)  = F_m%g(pairOne)  + F_m%g(pairTwo)
  F_m%rh(pairOne) = F_m%rh(pairOne) + F_m%rh(pairTwo)
  F_m%g(pairTwo)  = 1.0d0
  F_m%rh(pairTwo) = 0.0d0

  elemconKey%ints(1) = dst
  elemconKey%ints(2) = src
  CALL elemcon%get(elemconKey, pairThree)

  F_m%g(pairThree) = -1.0d0
enddo

deallocate(elemconKey%ints)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine FemApplyPeriodicity
