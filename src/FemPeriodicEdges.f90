subroutine FemPeriodicEdges(nodePairingFirst, nodePairingSecond, nodePairingSecondInverse, elemcon)
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod, only: F_m
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: nodePairingFirst, nodePairingSecond, nodePairingSecondInverse
type(fhash_type_iterator__ints_double)       :: nodePairingSecondIt
type(ints_type)                              :: nodePairingSecondKey, nodePairingSecondInverseKey
type(ints_type)                              :: destBothKey
integer                                      :: nodePairingSecondValue

type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemconKey

integer :: source2, dest2, dest1
integer :: destBoth
integer :: ii, mm, nn

logical :: success
!------------------------------------------------------------------------------------------------------!
allocate(elemconKey%ints(2))
allocate(destBothKey%ints(1))
allocate(nodePairingSecondInverseKey%ints(1))

call nodePairingSecondIt%begin(nodePairingSecond)
do ii = 1, nodePairingSecond%key_count()
  call nodePairingSecondIt%next(nodePairingSecondKey, nodePairingSecondValue)
  source2 = nodePairingSecondKey%ints(1)  ! source2 = 7 = source1
  dest2   = nodePairingSecondValue        ! dest2   = 1

  destBothKey%ints(1) = dest2
  call nodePairingFirst%get(destBothKey, destBoth, success) ! destBoth = 38

  if (success) then
    ! if success, we want to find the dest1 which has source1 = source2
    ! we can find it from the inverse hash as follows:
    ! dest1   = inverse_hash_yy(destBoth)        ! 23 = inverse_hash_yy(38) = hash_xx(7)
    ! source1 = source2 = inverse_hash_xx(dest1) ! 7

    nodePairingSecondInverseKey%ints(1) = destBoth
    call nodePairingSecondInverse%get(nodePairingSecondInverseKey, dest1)

    elemconKey%ints(1) = dest1
    elemconKey%ints(2) = source2
    call elemcon%get(elemconKey, mm) ! id of pair (23,7)

    elemconKey%ints(1) = destBoth
    elemconKey%ints(2) = dest2
    call elemcon%get(elemconKey, nn) ! id of pair (38,1)

    F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
    F_m%g(nn) = 0.0d0
  else
    cycle
  endif
enddo

deallocate(destBothKey%ints)
deallocate(elemconKey%ints)
deallocate(nodePairingSecondInverseKey%ints)

return
!------------------------------------------------------------------------------------------------------!
end subroutine FemPeriodicEdges
