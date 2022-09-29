subroutine FemPeriodicEdges(elemcon)
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod, only: F_m
use geometry_mod, only: edgeNodeOne, edgeNodeTwo, edgeNodeThree, edgeNodeFour, numEdgePeriodicPairs
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemconKey

integer :: ii, mm, nn
!------------------------------------------------------------------------------------------------------!
allocate(elemconKey%ints(2))

do ii = 1, numEdgePeriodicPairs
  elemconKey%ints(1) = edgeNodeOne(ii)
  elemconKey%ints(2) = edgeNodeTwo(ii)
  call elemcon%get(elemconKey, mm) ! id of pair (23,7)

  elemconKey%ints(1) = edgeNodeThree(ii)
  elemconKey%ints(2) = edgeNodeFour(ii)
  call elemcon%get(elemconKey, nn) ! id of pair (38,1)

  F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
  F_m%g(nn) = 0.0d0
enddo

deallocate(elemconKey%ints)

return
!------------------------------------------------------------------------------------------------------!
end subroutine FemPeriodicEdges
