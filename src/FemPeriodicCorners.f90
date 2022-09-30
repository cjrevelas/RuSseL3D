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

integer :: ii, mm, nn
!------------------------------------------------------------------------------------------------------!
allocate(elemconKey%ints(2))

do ii = 1, numCornerPeriodicPairsYY
  elemconKey%ints(1) = cornerNodeOneYY(ii) ! 9  or 48
  elemconKey%ints(2) = cornerNodeTwoYY(ii) ! 21 or 31
  call elemcon%get(elemconKey, mm)

  elemconKey%ints(1) = cornerNodeThreeYY(ii) ! 1 or 38 = destYY
  elemconKey%ints(2) = cornerNodeFourYY(ii)  ! 7 or 23
  call elemcon%get(elemconKey, nn)

  F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
  F_m%g(nn) = 0.0d0
enddo

do ii = 1, numCornerPeriodicPairsXX
  elemconKey%ints(1) = cornerNodeOneXX(1)  ! 31
  elemconKey%ints(2) = cornerNodeTwoXX(1)  ! 21
  call elemcon%get(elemconKey, mm)

  elemconKey%ints(1) = cornerNodeThreeXX(1) ! 23
  elemconKey%ints(2) = cornerNodeFourXX(1)  ! 7
  call elemcon%get(elemconKey, nn)

  F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
  F_m%g(nn) = 0.0d0
enddo

deallocate(elemconKey%ints)

return
!------------------------------------------------------------------------------------------------------!
end subroutine FemPeriodicCorners
