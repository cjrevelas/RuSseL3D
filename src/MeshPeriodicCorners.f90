subroutine MeshPeriodicCorners()
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use geometry_mod, only: nodePairingZZhash,                   &
                        nodePairingXXhashInverse,            &
                        nodePairingYYhashInverse,            &
                        cornerNodeOneYY, cornerNodeTwoYY,    &
                        cornerNodeThreeYY, cornerNodeFourYY, &
                        cornerNodeOneXX, cornerNodeTwoXX,    &
                        cornerNodeThreeXX, cornerNodeFourXX, &
                        destTriple, numCornerPeriodicPairsXX,&
                        numCornerPeriodicPairsYY
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
type(fhash_type_iterator__ints_double) :: nodePairingZZit
type(ints_type)                        :: nodePairingZZhashKey
type(ints_type)                        :: nodePairingXXhashInverseKey, nodePairingYYhashInverseKey
integer                                :: nodePairingZZvalue

logical :: success, successAux

integer :: sourceXX
integer :: sourceYY
integer :: sourceZZ, destZZ
integer :: sourceAux = 0
integer :: ii, kk
!------------------------------------------------------------------------------------------------------!
allocate(nodePairingYYhashInverseKey%ints(1))
allocate(nodePairingXXhashInverseKey%ints(1))

call nodePairingZZit%begin(nodePairingZZhash)

kk = 0
do ii = 1, nodePairingZZhash%key_count()
  call NodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)

  sourceZZ = nodePairingZZhashKey%ints(1)
  destZZ   = nodePairingZZvalue

  nodePairingYYhashInverseKey%ints(1) = destZZ
  call nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, sourceYY, success)

  if (success) kk = kk + 1
enddo

numCornerPeriodicPairsYY = kk
!write(123,*) "corner yy zz pairs:" , numCornerPeriodicPairsYY

allocate(cornerNodeOneYY(numCornerPeriodicPairsYY))
allocate(cornerNodeTwoYY(numCornerPeriodicPairsYY))
allocate(cornerNodeThreeYY(numCornerPeriodicPairsYY))
allocate(cornerNodeFourYY(numCornerPeriodicPairsYY))

call nodePairingZZit%begin(nodePairingZZhash)

kk = 0
do ii = 1, nodePairingZZhash%key_count()
  call nodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)

  sourceZZ = nodePairingZZhashKey%ints(1)
  destZZ   = nodePairingZZvalue

  nodePairingYYhashInverseKey%ints(1) = destZZ
  call nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, sourceYY, success)

  if (success) then
    kk = kk + 1
    nodePairingYYhashInverseKey%ints(1) = sourceZZ
    call nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, sourceAux)

    cornerNodeOneYY(kk) = sourceZZ
    cornerNodeTwoYY(kk) = sourceAux

    cornerNodeThreeYY(kk) = destZZ
    cornerNodeFourYY(kk)  = sourceYY

    nodePairingXXhashInverseKey%ints(1) = destZZ
    call nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceXX, successAux)
    if (successAux) destTriple = destZZ
  endif
enddo

call nodePairingZZit%begin(nodePairingZZhash)

kk = 0
do ii = 1, nodePairingZZhash%key_count()
  call nodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)
  sourceZZ = nodePairingZZhashKey%ints(1)
  destZZ   = nodePairingZZvalue

  nodePairingXXhashInverseKey%ints(1) = destZZ
  call nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceXX, success)

  if (success) kk = kk + 1
enddo

numCornerPeriodicPairsXX = kk
!write(123,*) "corner xx zz pairs:" , numCornerPeriodicPairsXX

allocate(cornerNodeOneXX(numCornerPeriodicPairsXX))
allocate(cornerNodeTwoXX(numCornerPeriodicPairsXX))
allocate(cornerNodeThreeXX(numCornerPeriodicPairsXX))
allocate(cornerNodeFourXX(numCornerPeriodicPairsXX))

call nodePairingZZit%begin(nodePairingZZhash)

kk = 0
do ii = 1, nodePairingZZhash%key_count()
  call nodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)

  sourceZZ = nodePairingZZhashKey%ints(1)
  destZZ   = nodePairingZZvalue

  nodePairingXXhashInverseKey%ints(1) = destZZ
  call nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceXX, success)

  if (success.and.(destZZ.ne.destTriple)) then
    kk = kk + 1

    nodePairingXXhashInverseKey%ints(1) = sourceZZ
    call nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceAux)

    cornerNodeOneXX(kk) = sourceZZ
    cornerNodeTwoXX(kk) = sourceAux

    cornerNodeThreeXX(kk) = destZZ
    cornerNodeFourXX(kk)  = sourceXX
  endif
enddo

deallocate(nodePairingYYhashInverseKey%ints)
deallocate(nodePairingXXhashInverseKey%ints)

return
!------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicCorners
