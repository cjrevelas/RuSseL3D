subroutine FemPeriodicEdges(elemcon)
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod, only: F_m
use parser_vars_mod, only: periodicAxisId
use geometry_mod,    only: edgeNodeOneXY, edgeNodeTwoXY, edgeNodeThreeXY, edgeNodeFourXY, numEdgePeriodicPairsXY, &
                           edgeNodeOneXZ, edgeNodeTwoXZ, edgeNodeThreeXZ, edgeNodeFourXZ, numEdgePeriodicPairsXZ, &
                           edgeNodeOneYZ, edgeNodeTwoYZ, edgeNodeThreeYZ, edgeNodeFourYZ, numEdgePeriodicPairsYZ
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemconKey

integer :: edge, pairOne, pairTwo
!------------------------------------------------------------------------------------------------------!
allocate(elemconKey%ints(2))

! Apply correction to periodic node pairs belonging to XY edges
if (periodicAxisId(1).and.periodicAxisId(2)) then
  do edge = 1, numEdgePeriodicPairsXY
    elemconKey%ints(1) = edgeNodeOneXY(edge)
    elemconKey%ints(2) = edgeNodeTwoXY(edge)
    call elemcon%get(elemconKey, pairOne) ! id of pair (23,7)

    elemconKey%ints(1) = edgeNodeThreeXY(edge)
    elemconKey%ints(2) = edgeNodeFourXY(edge)
    call elemcon%get(elemconKey, pairTwo) ! id of pair (38,1)

    F_m%g(pairOne) = F_m%g(pairOne) + F_m%g(pairTwo)
    F_m%g(pairTwo) = 0.0d0
  enddo
endif

! Apply correction to periodic node pairs belonging to XZ edges
if (periodicAxisId(1).and.periodicAxisId(3)) then
  do edge = 1, numEdgePeriodicPairsXZ
    elemconKey%ints(1) = edgeNodeOneXZ(edge)
    elemconKey%ints(2) = edgeNodeTwoXZ(edge)
    call elemcon%get(elemconKey, pairOne)

    elemconKey%ints(1) = edgeNodeThreeXZ(edge)
    elemconKey%ints(2) = edgeNodeFourXZ(edge)
    call elemcon%get(elemconKey, pairTwo)

    F_m%g(pairOne) = F_m%g(pairOne) + F_m%g(pairTwo)
    F_m%g(pairTwo) = 0.0d0
  enddo
endif

! Apply correction to periodic node pairs belonging to YZ edges
if (periodicAxisId(2).and.periodicAxisId(3)) then
  do edge = 1, numEdgePeriodicPairsYZ
    elemconKey%ints(1) = edgeNodeOneYZ(edge)
    elemconKey%ints(2) = edgeNodeTwoYZ(edge)
    call elemcon%get(elemconKey, pairOne)

    elemconKey%ints(1) = edgeNodeThreeYZ(edge)
    elemconKey%ints(2) = edgeNodeFourYZ(edge)
    call elemcon%get(elemconKey, pairTwo)

    F_m%g(pairOne) = F_m%g(pairOne) + F_m%g(pairTwo)
    F_m%g(pairTwo) = 0.0d0
  enddo
endif

deallocate(elemconKey%ints)

return
!------------------------------------------------------------------------------------------------------!
end subroutine FemPeriodicEdges
