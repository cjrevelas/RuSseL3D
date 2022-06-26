subroutine FemPeriodicEdges(nodePairingFirst, nodePairingSecond, elemcon)
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod, only: F_m
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: nodePairingFirst, nodePairingSecond
type(fhash_type_iterator__ints_double)       :: nodePairingFirstIt, nodePairingSecondIt
type(ints_type)                              :: nodePairingFirstKey, nodePairingSecondKey
type(ints_type)                              :: destBothKey
integer                                      :: nodePairingFirstValue, nodePairingSecondValue

type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemcon_key

integer :: source1, dest1
integer :: source2, dest2
integer :: destBoth
integer :: ii, jj, mm, nn

logical :: success
!------------------------------------------------------------------------------------------------------!
call nodePairingSecondIt%begin(nodePairingSecond)

allocate(elemcon_key%ints(2))

allocate(destBothKey%ints(1))

do ii = 1, nodePairingSecond%key_count()
  call nodePairingSecondIt%next(nodePairingSecondKey, nodePairingSecondValue)
  source2 = nodePairingSecondKey%ints(1)  ! source2 = 7
  dest2   = nodePairingSecondValue        ! dest2   = 1

  destBothKey%ints(1) = dest2
  call nodePairingFirst%get(destBothKey, destBoth, success) ! destBoth = 38

  call nodePairingFirstIt%begin(nodePairingFirst)
  do jj = 1, nodePairingFirst%key_count()
    call nodePairingFirstIt%next(nodePairingFirstKey, nodePairingFirstValue)
    source1 = nodePairingFirstKey%ints(1) ! source1 = 7
    dest1   = nodePairingFirstValue       ! dest1   = 23

    if (source2.eq.source1) then
      elemcon_key%ints(1) = dest1
      elemcon_key%ints(2) = source1
      call elemcon%get(elemcon_key, mm) ! id of pair (23,7)

      elemcon_key%ints(1) = destBoth
      elemcon_key%ints(2) = dest2
      call elemcon%get(elemcon_key, nn) ! id of pair (38,1)

      F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
      F_m%g(nn) = 0.0d0
    endif
  enddo
enddo

deallocate(destBothKey%ints)
deallocate(elemcon_key%ints)
return
!------------------------------------------------------------------------------------------------------!
end subroutine FemPeriodicEdges
