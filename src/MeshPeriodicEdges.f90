subroutine MeshPeriodicEdges(nodePairingFirst, nodePairingSecond, nodePairingSecondInverse)
!----------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use geometry_mod, only: edgeNodeOne, edgeNodeTwo, edgeNodeThree, edgeNodeFour, &
                        numEdgePeriodicPairs
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!

type(fhash_type__ints_double), intent(inout) :: nodePairingFirst, nodePairingSecond
type(fhash_type__ints_double), intent(inout) :: nodePairingSecondInverse
type(fhash_type_iterator__ints_double)       :: nodePairingSecondIt
type(ints_type)                              :: nodePairingSecondKey
type(ints_type)                              :: nodePairingSecondInverseKey
type(ints_type)                              :: destBothKey
integer                                      :: nodePairingSecondValue

integer :: dest1, source2, destBoth, dest2, ii, kk

logical :: success
!----------------------------------------------------------------------------------------------------------------------------!
allocate(destBothKey%ints(1))
allocate(nodePairingSecondInverseKey%ints(1))

kk = 0 ! counter of periodic pairs belonging to edges and thus need correction later
call nodePairingSecondIt%begin(nodePairingSecond)
do ii = 1, nodePairingSecond%key_count()
  call nodePairingSecondIt%next(nodePairingSecondKey, nodePairingSecondValue)

  source2 = nodePairingSecondKey%ints(1)  ! source2 = 7 = source1
  dest2   = nodePairingSecondValue        ! dest2   = 1

  destBothKey%ints(1) = dest2
  call nodePairingFirst%get(destBothKey, destBoth, success)  ! destBoth = 38

  if (success) kk = kk + 1
enddo

numEdgePeriodicPairs = kk

allocate(edgeNodeOne(numEdgePeriodicPairs))
allocate(edgeNodeTwo(numEdgePeriodicPairs))
allocate(edgeNodeThree(numEdgePeriodicPairs))
allocate(edgeNodeFour(numEdgePeriodicPairs))

kk = 0
call nodePairingSecondIt%begin(nodePairingSecond)
do ii = 1, nodePairingSecond%key_count()
  call nodePairingSecondIt%next(nodePairingSecondKey, nodePairingSecondValue)

  source2 = nodePairingSecondKey%ints(1)  ! source2 = 7 = source1
  dest2   = nodePairingSecondValue        ! dest2   = 1

  destBothKey%ints(1) = dest2
  call nodePairingFirst%get(destBothKey, destBoth, success)  ! destBoth = 38

  if (success) then
    kk = kk + 1

    ! if success, we want to find the dest1 which has source1 = source2
    ! we can find it from the inverse hash as follows:
    ! dest1   = inverse_hash_yy(destBoth)        ! 23 = inverse_hash_yy(38) = hash_xx(7)
    ! source1 = source2 = inverse_hash_xx(dest1) ! 7

    nodePairingSecondInverseKey%ints(1) = destBoth
    call nodePairingSecondInverse%get(nodePairingSecondInverseKey, dest1)

    edgeNodeOne(kk) = dest1
    edgeNodeTwo(kk) = source2

    edgeNodeThree(kk) = destBoth
    edgeNodeFour(kk)  = dest2
  else
    cycle
  endif
enddo

deallocate(destBothKey%ints)
deallocate(nodePairingSecondInverseKey%ints)

return
!----------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicEdges
