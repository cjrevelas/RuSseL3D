subroutine FemPeriodicCorners(elemcon)
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod,      only: F_m
use geometry_mod, only: cornerNodeOneXX, cornerNodeTwoXX, cornerNodeThreeXX, cornerNodeFourXX, &
                        cornerNodeOneYY, cornerNodeTwoYY, cornerNodeThreeYY, cornerNodeFourYY, &
                        numCornerPeriodicPairsXX, numCornerPeriodicPairsYY
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemconKey

integer :: corner, pairOne, pairTwo
!------------------------------------------------------------------------------------------------------!
allocate(elemconKey%ints(2))

do corner = 1, numCornerPeriodicPairsYY
  elemconKey%ints(1) = cornerNodeOneYY(corner) ! 9  or 48
  elemconKey%ints(2) = cornerNodeTwoYY(corner) ! 21 or 31
  CALL elemcon%get(elemconKey, pairOne)

  elemconKey%ints(1) = cornerNodeThreeYY(corner) ! 1 or 38 = dstYY
  elemconKey%ints(2) = cornerNodeFourYY(corner)  ! 7 or 23
  CALL elemcon%get(elemconKey, pairTwo)

  F_m%g(pairOne) = F_m%g(pairOne) + F_m%g(pairTwo)
  F_m%g(pairTwo) = 0.0d0
enddo

do corner = 1, numCornerPeriodicPairsXX
  elemconKey%ints(1) = cornerNodeOneXX(1)  ! 31
  elemconKey%ints(2) = cornerNodeTwoXX(1)  ! 21
  CALL elemcon%get(elemconKey, pairOne)

  elemconKey%ints(1) = cornerNodeThreeXX(1) ! 23
  elemconKey%ints(2) = cornerNodeFourXX(1)  ! 7
  CALL elemcon%get(elemconKey, pairTwo)

  F_m%g(pairOne) = F_m%g(pairOne) + F_m%g(pairTwo)
  F_m%g(pairTwo) = 0.0d0
enddo

deallocate(elemconKey%ints)

return
!------------------------------------------------------------------------------------------------------!
end subroutine FemPeriodicCorners
