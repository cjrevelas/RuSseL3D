!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshAppendPeriodicPairs(elemcon, startingPair, endingPair, nodePairingHash)
!----------------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod,      only: F_m
use geometry_mod, only: nodePairId
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemconKey
integer                                      :: elemconValue

integer, intent(in) :: startingPair, endingPair

type(fhash_type__ints_double), intent(inout) :: nodePairingHash
type(fhash_type_iterator__ints_double)       :: nodePairingIt
type(ints_type)                              :: nodePairingKey
integer                                      :: nodePairingValue

integer :: source, dest
integer :: kk

logical :: success
!----------------------------------------------------------------------------------------------------------------------------------!
allocate(elemconKey%ints(2))

call nodePairingIt%begin(nodePairingHash)

do kk = startingPair, endingPair
  call nodePairingIt%next(nodePairingKey, nodePairingValue)

  source = nodePairingKey%ints(1)
  dest   = nodePairingValue

  ! Append pair
  F_m%row(kk) = dest
  F_m%col(kk) = source

  elemconKey%ints(1) = dest
  elemconKey%ints(2) = source

  call elemcon%get(elemconKey, elemconValue, success)

  if (success) then
    nodePairId(kk) = elemconValue
  else
    call elemcon%set(elemconKey, kk)
    nodePairId(kk) = kk
  endif
enddo

deallocate(elemconKey%ints)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshAppendPeriodicPairs
