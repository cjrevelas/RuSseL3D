subroutine FemPeriodicCorners(elemcon)
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use kcw_mod,      only: F_m
use geometry_mod, only: nodePairingZZhash,        &
                        nodePairingXXhashInverse, &
                        nodePairingYYhashInverse
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
type(fhash_type_iterator__ints_double) :: nodePairingZZit
type(ints_type)                        :: nodePairingZZhashKey
type(ints_type)                        :: nodePairingXXhashInverseKey, nodePairingYYhashInverseKey
integer                                :: nodePairingZZvalue

type(fhash_type__ints_double), intent(inout) :: elemcon
type(ints_type)                              :: elemconKey

logical :: success, successAux

integer :: sourceXX
integer :: sourceYY
integer :: sourceZZ, destZZ
integer :: sourceAux=0, destTriple=0
integer :: ii, mm, nn
!------------------------------------------------------------------------------------------------------!
allocate(elemconKey%ints(2))
allocate(nodePairingYYhashInverseKey%ints(1))
allocate(nodePairingXXhashInverseKey%ints(1))

call nodePairingZZit%begin(nodePairingZZhash)

do ii = 1, nodePairingZZhash%key_count()
  call nodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)
  sourceZZ = nodePairingZZhashKey%ints(1) ! sourceZZ = 9 or 48
  destZZ   = nodePairingZZvalue           ! destZZ   = 1 or 38

  nodePairingYYhashInverseKey%ints(1) = destZZ
  call nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, sourceYY, success) ! sourceAux = 7 or 23

  if (success) then
    nodePairingYYhashInverseKey%ints(1) = sourceZZ
    call nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, sourceAux) ! 21 or 31

    elemconKey%ints(1) = sourceZZ  ! 9  or 48
    elemconKey%ints(2) = sourceAux ! 21 or 31
    call elemcon%get(elemconKey, mm)

    elemconKey%ints(1) = destZZ   ! 1 or 38 = destYY
    elemconKey%ints(2) = sourceYY ! 7 or 23
    call elemcon%get(elemconKey, nn)

    F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
    F_m%g(nn) = 0.0d0

    nodePairingXXhashInverseKey%ints(1) = destZZ ! = destYY = 1 or 38
    call nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceXX, successAux)
    if (successAux) destTriple = destZZ ! 38
  endif
enddo

call nodePairingZZit%begin(nodePairingZZhash)
do ii = 1, nodePairingZZhash%key_count()
  call nodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)
  sourceZZ = nodePairingZZhashKey%ints(1) ! sourceZZ = 31
  destZZ   = nodePairingZZvalue           ! destZZ   = 23

  nodePairingXXhashInverseKey%ints(1) = destZZ
  call nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceXX, success) ! sourceXX = 7

  if (success.and.(destZZ.ne.destTriple)) then
    nodePairingXXhashInverseKey%ints(1) = sourceZZ
    call nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceAux) ! sourceAux = 21

    elemconKey%ints(1) = sourceZZ   ! 31
    elemconKey%ints(2) = sourceAux  ! 21
    call elemcon%get(elemconKey, mm)

    elemconKey%ints(1) = destZZ     ! 23
    elemconKey%ints(2) = sourceXX   ! 7
    call elemcon%get(elemconKey, nn)

    F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
    F_m%g(nn) = 0.0d0
  endif
enddo

deallocate(elemconKey%ints)
deallocate(nodePairingYYhashInverseKey%ints)
deallocate(nodePairingXXhashInverseKey%ints)

return
!------------------------------------------------------------------------------------------------------!
end subroutine FemPeriodicCorners
