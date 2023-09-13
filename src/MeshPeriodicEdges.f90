subroutine MeshPeriodicEdges()
!----------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: periodicAxisId
use geometry_mod,    only: edgeNodeOneXY, edgeNodeTwoXY, edgeNodeThreeXY, edgeNodeFourXY, numEdgePeriodicPairsXY, &
                           edgeNodeOneXZ, edgeNodeTwoXZ, edgeNodeThreeXZ, edgeNodeFourXZ, numEdgePeriodicPairsXZ, &
                           edgeNodeOneYZ, edgeNodeTwoYZ, edgeNodeThreeYZ, edgeNodeFourYZ, numEdgePeriodicPairsYZ, &
                           nodePairingYYhash, nodePairingYYhashInverse,                                           &
                           nodePairingZZhash, nodePairingZZhashInverse,                                           &
                           nodePairingXXhash
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
type(fhash_type_iterator__ints_double) :: nodePairingYYit, nodePairingZZit
type(ints_type)                        :: nodePairingYYkey, nodePairingZZkey
type(ints_type)                        :: nodePairingYYinverseKey, nodePairingZZinverseKey
type(ints_type)                        :: destBothKey
integer                                :: nodePairingYYvalue, nodePairingZZvalue

integer :: dest1, source2, destBoth, dest2, keyIndex, pair, numPairs

logical :: success
!----------------------------------------------------------------------------------------------------------------------------!
allocate(destBothKey%ints(1))
allocate(nodePairingYYinverseKey%ints(1))
allocate(nodePairingZZinverseKey%ints(1))

! Find number of periodic node pairs between XY edges
if (periodicAxisId(1).and.periodicAxisId(2)) then
  numPairs = 0

  call nodePairingYYit%begin(nodePairingYYhash)
  do keyIndex = 1, nodePairingYYhash%key_count()
    call nodePairingYYit%next(nodePairingYYkey, nodePairingYYvalue)

    source2 = nodePairingYYkey%ints(1) ! source2 = 7 = source1
    dest2   = nodePairingYYvalue       ! dest2   = 1

    destBothKey%ints(1) = dest2
    call nodePairingXXhash%get(destBothKey, destBoth, success)  ! destBoth = 38

    if (success) numPairs = numPairs + 1
  enddo

  numEdgePeriodicPairsXY = numPairs

  allocate(edgeNodeOneXY(numEdgePeriodicPairsXY))
  allocate(edgeNodeTwoXY(numEdgePeriodicPairsXY))
  allocate(edgeNodeThreeXY(numEdgePeriodicPairsXY))
  allocate(edgeNodeFourXY(numEdgePeriodicPairsXY))
endif

! Find number of periodic node pairs between XZ edges
if (periodicAxisId(1).and.periodicAxisId(3)) then
  numPairs = 0

  call nodePairingZZit%begin(nodePairingZZhash)
  do keyIndex = 1, nodePairingZZhash%key_count()
    call nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)

    source2 = nodePairingZZkey%ints(1) ! source2 = 7 = source1
    dest2   = nodePairingZZvalue       ! dest2   = 1

    destBothKey%ints(1) = dest2
    call nodePairingXXhash%get(destBothKey, destBoth, success)  ! destBoth = 38

    if (success) numPairs = numPairs + 1
  enddo

  numEdgePeriodicPairsXZ = numPairs

  allocate(edgeNodeOneXZ(numEdgePeriodicPairsXZ))
  allocate(edgeNodeTwoXZ(numEdgePeriodicPairsXZ))
  allocate(edgeNodeThreeXZ(numEdgePeriodicPairsXZ))
  allocate(edgeNodeFourXZ(numEdgePeriodicPairsXZ))
endif

! Find number of periodic node pairs between YZ edges
if (periodicAxisId(2).and.periodicAxisId(3)) then
  numPairs = 0

  call nodePairingZZit%begin(nodePairingZZhash)
  do keyIndex = 1, nodePairingZZhash%key_count()
    call nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)

    source2 = nodePairingZZkey%ints(1) ! source2 = 7 = source1
    dest2   = nodePairingZZvalue       ! dest2   = 1

    destBothKey%ints(1) = dest2
    call nodePairingYYhash%get(destBothKey, destBoth, success)  ! destBoth = 38

    if (success) numPairs = numPairs + 1
  enddo

  numEdgePeriodicPairsYZ = numPairs

  allocate(edgeNodeOneYZ(numEdgePeriodicPairsYZ))
  allocate(edgeNodeTwoYZ(numEdgePeriodicPairsYZ))
  allocate(edgeNodeThreeYZ(numEdgePeriodicPairsYZ))
  allocate(edgeNodeFourYZ(numEdgePeriodicPairsYZ))
endif

! Perform XY node pair association
if (periodicAxisId(1).and.periodicAxisId(2)) then
  pair = 0

  call nodePairingYYit%begin(nodePairingYYhash)
  do keyIndex = 1, nodePairingYYhash%key_count()
    call nodePairingYYit%next(nodePairingYYkey, nodePairingYYvalue)

    source2 = nodePairingYYkey%ints(1)  ! source2 = 7 = source1
    dest2   = nodePairingYYvalue        ! dest2   = 1

    destBothKey%ints(1) = dest2
    call nodePairingXXhash%get(destBothKey, destBoth, success)  ! destBoth = 38

    if (success) then
      pair = pair + 1

      ! if success, we want to find the dest1 which has source1 = source2
      ! we can find it from the inverse hash as follows:
      ! dest1   = inverse_hash_yy(destBoth)        ! 23 = inverse_hash_yy(38) = hash_xx(7)
      ! source1 = source2 = inverse_hash_xx(dest1) ! 7

      nodePairingYYinverseKey%ints(1) = destBoth
      call nodePairingYYhashInverse%get(nodePairingYYinverseKey, dest1)

      edgeNodeOneXY(pair) = dest1
      edgeNodeTwoXY(pair) = source2

      edgeNodeThreeXY(pair) = destBoth
      edgeNodeFourXY(pair)  = dest2
    else
      cycle
    endif
  enddo
endif


! Perform XZ node pair association
if (periodicAxisId(1).and.periodicAxisId(3)) then
  pair = 0

  call nodePairingZZit%begin(nodePairingZZhash)
  do keyIndex = 1, nodePairingZZhash%key_count()
    call nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)

    source2 = nodePairingZZkey%ints(1)  ! source2 = 7 = source1
    dest2   = nodePairingZZvalue        ! dest2   = 1

    destBothKey%ints(1) = dest2
    call nodePairingXXhash%get(destBothKey, destBoth, success)  ! destBoth = 38

    if (success) then
      pair = pair + 1

      nodePairingZZinverseKey%ints(1) = destBoth
      call nodePairingZZhashInverse%get(nodePairingZZinverseKey, dest1)

      edgeNodeOneXZ(pair) = dest1
      edgeNodeTwoXZ(pair) = source2

      edgeNodeThreeXZ(pair) = destBoth
      edgeNodeFourXZ(pair)  = dest2
    else
      cycle
    endif
  enddo
endif

! Perform YZ node pair association
if (periodicAxisId(2).and.periodicAxisId(3)) then
  pair = 0

  call nodePairingZZit%begin(nodePairingZZhash)
  do keyIndex = 1, nodePairingZZhash%key_count()
    call nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)

    source2 = nodePairingZZkey%ints(1)  ! source2 = 7 = source1
    dest2   = nodePairingZZvalue        ! dest2   = 1

    destBothKey%ints(1) = dest2
    call nodePairingYYhash%get(destBothKey, destBoth, success)  ! destBoth = 38

    if (success) then
      pair = pair + 1

      nodePairingZZinverseKey%ints(1) = destBoth
      call nodePairingZZhashInverse%get(nodePairingZZinverseKey, dest1)

      edgeNodeOneYZ(pair) = dest1
      edgeNodeTwoYZ(pair) = source2

      edgeNodeThreeYZ(pair) = destBoth
      edgeNodeFourYZ(pair)  = dest2
    else
      cycle
    endif
  enddo
endif

deallocate(destBothKey%ints)
deallocate(nodePairingYYinverseKey%ints)
deallocate(nodePairingZZinverseKey%ints)

return
!----------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicEdges
