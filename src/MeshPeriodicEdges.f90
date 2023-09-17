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
type(ints_type)                        :: dstBothKey
integer                                :: nodePairingYYvalue, nodePairingZZvalue

integer :: dst1, src2, dstBoth, dst2, keyIndex, pair, numPairs

logical :: success
!----------------------------------------------------------------------------------------------------------------------------!
allocate(dstBothKey%ints(1))
allocate(nodePairingYYinverseKey%ints(1))
allocate(nodePairingZZinverseKey%ints(1))

! Find number of periodic node pairs between XY edges
if (periodicAxisId(1).and.periodicAxisId(2)) then
  numPairs = 0

  CALL nodePairingYYit%begin(nodePairingYYhash)
  do keyIndex = 1, nodePairingYYhash%key_count()
    CALL nodePairingYYit%next(nodePairingYYkey, nodePairingYYvalue)

    src2 = nodePairingYYkey%ints(1) ! src2 = 7 = src1
    dst2 = nodePairingYYvalue       ! dst2 = 1

    dstBothKey%ints(1) = dst2
    CALL nodePairingXXhash%get(dstBothKey, dstBoth, success) ! dstBoth = 38

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

  CALL nodePairingZZit%begin(nodePairingZZhash)
  do keyIndex = 1, nodePairingZZhash%key_count()
    CALL nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)

    src2 = nodePairingZZkey%ints(1) ! src2 = 7 = src1
    dst2 = nodePairingZZvalue       ! dst2 = 1

    dstBothKey%ints(1) = dst2
    CALL nodePairingXXhash%get(dstBothKey, dstBoth, success) ! dstBoth = 38

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

  CALL nodePairingZZit%begin(nodePairingZZhash)
  do keyIndex = 1, nodePairingZZhash%key_count()
    CALL nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)

    src2 = nodePairingZZkey%ints(1) ! src2 = 7 = src1
    dst2 = nodePairingZZvalue       ! dst2 = 1

    dstBothKey%ints(1) = dst2
    CALL nodePairingYYhash%get(dstBothKey, dstBoth, success) ! dstBoth = 38

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

  CALL nodePairingYYit%begin(nodePairingYYhash)
  do keyIndex = 1, nodePairingYYhash%key_count()
    CALL nodePairingYYit%next(nodePairingYYkey, nodePairingYYvalue)

    src2 = nodePairingYYkey%ints(1)  ! src2 = 7 = src1
    dst2 = nodePairingYYvalue        ! dst2 = 1

    dstBothKey%ints(1) = dst2
    CALL nodePairingXXhash%get(dstBothKey, dstBoth, success)  ! dstBoth = 38

    if (success) then
      pair = pair + 1

      ! if success, we want to find the dst1 which has src1 = src2
      ! we can find it from the inverse hash as follows:
      ! dst1 = inverseHashYY(dstBoth)     ! 23 = inverseHashYY(38) = hashXX(7)
      ! src1 = src2 = inverseHashXX(dst1) ! 7

      nodePairingYYinverseKey%ints(1) = dstBoth
      CALL nodePairingYYhashInverse%get(nodePairingYYinverseKey, dst1)

      edgeNodeOneXY(pair) = dst1
      edgeNodeTwoXY(pair) = src2

      edgeNodeThreeXY(pair) = dstBoth
      edgeNodeFourXY(pair)  = dst2
    else
      cycle
    endif
  enddo
endif

! Perform XZ node pair association
if (periodicAxisId(1).and.periodicAxisId(3)) then
  pair = 0

  CALL nodePairingZZit%begin(nodePairingZZhash)
  do keyIndex = 1, nodePairingZZhash%key_count()
    CALL nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)

    src2 = nodePairingZZkey%ints(1)  ! src2 = 7 = src1
    dst2 = nodePairingZZvalue        ! dst2 = 1

    dstBothKey%ints(1) = dst2
    CALL nodePairingXXhash%get(dstBothKey, dstBoth, success) ! dstBoth = 38

    if (success) then
      pair = pair + 1

      nodePairingZZinverseKey%ints(1) = dstBoth
      CALL nodePairingZZhashInverse%get(nodePairingZZinverseKey, dst1)

      edgeNodeOneXZ(pair) = dst1
      edgeNodeTwoXZ(pair) = src2

      edgeNodeThreeXZ(pair) = dstBoth
      edgeNodeFourXZ(pair)  = dst2
    else
      cycle
    endif
  enddo
endif

! Perform YZ node pair association
if (periodicAxisId(2).and.periodicAxisId(3)) then
  pair = 0

  CALL nodePairingZZit%begin(nodePairingZZhash)
  do keyIndex = 1, nodePairingZZhash%key_count()
    CALL nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)

    src2 = nodePairingZZkey%ints(1) ! src2 = 7 = src1
    dst2 = nodePairingZZvalue       ! dst2 = 1

    dstBothKey%ints(1) = dst2
    CALL nodePairingYYhash%get(dstBothKey, dstBoth, success) ! dstBoth = 38

    if (success) then
      pair = pair + 1

      nodePairingZZinverseKey%ints(1) = dstBoth
      CALL nodePairingZZhashInverse%get(nodePairingZZinverseKey, dst1)

      edgeNodeOneYZ(pair) = dst1
      edgeNodeTwoYZ(pair) = src2

      edgeNodeThreeYZ(pair) = dstBoth
      edgeNodeFourYZ(pair)  = dst2
    else
      cycle
    endif
  enddo
endif

deallocate(dstBothKey%ints)
deallocate(nodePairingYYinverseKey%ints)
deallocate(nodePairingZZinverseKey%ints)

return
!----------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicEdges
