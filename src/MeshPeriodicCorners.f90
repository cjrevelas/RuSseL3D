subroutine MeshPeriodicCorners()
!------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
!use constants_mod, only: tol
use geometry_mod, only: nodePairingZZhash,                    &
                        nodePairingXXhashInverse,             &
                        nodePairingYYhashInverse,             &
                        cornerNodeOneYY, cornerNodeTwoYY,     &
                        cornerNodeThreeYY, cornerNodeFourYY,  &
                        cornerNodeOneXX, cornerNodeTwoXX,     &
                        cornerNodeThreeXX, cornerNodeFourXX,  &
                        dstTriple, numCornerPeriodicPairsXX, &
                        numCornerPeriodicPairsYY!, nodeCoord,  &
                        !boxLow, boxHigh
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
type(fhash_type_iterator__ints_double) :: nodePairingZZit
type(ints_type)                        :: nodePairingZZhashKey
type(ints_type)                        :: nodePairingXXhashInverseKey, nodePairingYYhashInverseKey
integer                                :: nodePairingZZvalue

logical :: success, successAux

integer :: srcXX
integer :: srcYY
integer :: srcZZ, dstZZ
integer :: srcAux = 0
integer :: counter, keyIndex
!------------------------------------------------------------------------------------------------------!
allocate(nodePairingYYhashInverseKey%ints(1))
allocate(nodePairingXXhashInverseKey%ints(1))

CALL nodePairingZZit%begin(nodePairingZZhash)

counter = 0
do keyIndex = 1, nodePairingZZhash%key_count()
  CALL NodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)

  srcZZ = nodePairingZZhashKey%ints(1)
  dstZZ = nodePairingZZvalue

  nodePairingYYhashInverseKey%ints(1) = dstZZ
  CALL nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, srcYY, success)

  !if ((success).and.((abs(nodeCoord(1,dstZZ)-boxLow(1))<tol).or.(abs(nodeCoord(1,dstZZ)-boxHigh(1))<tol))) counter = counter + 1
  if (success) counter = counter + 1
enddo

numCornerPeriodicPairsYY = counter
!write(123,*) "corner yy zz pairs:", numCornerPeriodicPairsYY
!write(123,*) "-----------------------"

allocate(cornerNodeOneYY(numCornerPeriodicPairsYY))
allocate(cornerNodeTwoYY(numCornerPeriodicPairsYY))
allocate(cornerNodeThreeYY(numCornerPeriodicPairsYY))
allocate(cornerNodeFourYY(numCornerPeriodicPairsYY))

CALL nodePairingZZit%begin(nodePairingZZhash)

counter = 0
do keyIndex = 1, nodePairingZZhash%key_count()
  CALL nodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)

  srcZZ = nodePairingZZhashKey%ints(1)
  dstZZ = nodePairingZZvalue

  nodePairingYYhashInverseKey%ints(1) = dstZZ ! which is equal to dstZZ
  CALL nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, srcYY, success)

  !if ((success).and.((abs(nodeCoord(1,dstZZ)-boxLow(1))<tol).or.(abs(nodeCoord(1,dstZZ)-boxHigh(1))<tol))) then
  if (success) then
    counter = counter + 1

    nodePairingYYhashInverseKey%ints(1) = srcZZ
    CALL nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, srcAux)

    cornerNodeOneYY(counter) = srcZZ
    cornerNodeTwoYY(counter) = srcAux

    cornerNodeThreeYY(counter) = dstZZ
    cornerNodeFourYY(counter)  = srcYY

    !write(123,*) "srcZZ: ", srcZZ, " -> dstZZ: ", dstZZ,    " [", nodeCoord(1,srcZZ), ", ", nodeCoord(2,srcZZ), ", ", nodeCoord(3,srcZZ), "] -> [", nodeCoord(1,dstZZ),   ", ", nodeCoord(2,dstZZ),   ", ", nodeCoord(3,dstZZ),   "]"
    !write(123,*) "dstYY: ", dstZZ,   " -> srcYY: ", srcYY , " [", nodeCoord(1,dstZZ),   ", ", nodeCoord(2,dstZZ),   ", ", nodeCoord(3,dstZZ),   "] -> [", nodeCoord(1,srcYY), ", ", nodeCoord(2,srcYY), ", ", nodeCoord(3,srcYY), "]"
    !write(123,*) "-----------------------"

    nodePairingXXhashInverseKey%ints(1) = dstZZ
    CALL nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, srcXX, successAux)
    if (successAux) then
      dstTriple = dstZZ
      !write(123,*) "FOUND dstTriple: ", dstZZ
    endif
  endif
enddo

CALL nodePairingZZit%begin(nodePairingZZhash)

counter = 0
do keyIndex = 1, nodePairingZZhash%key_count()
  CALL nodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)
  srcZZ = nodePairingZZhashKey%ints(1)
  dstZZ = nodePairingZZvalue

  nodePairingXXhashInverseKey%ints(1) = dstZZ
  CALL nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, srcXX, success)

  !if ((success).and.((abs(nodeCoord(2,dstZZ)-boxLow(2))<tol).or.(abs(nodeCoord(2,dstZZ)-boxHigh(2))<tol))) counter = counter + 1
  if (success) counter = counter + 1
enddo

numCornerPeriodicPairsXX = counter
!write(456,*) "corner xx zz pairs:" , numCornerPeriodicPairsXX
!write(456,*) "-----------------------"

allocate(cornerNodeOneXX(numCornerPeriodicPairsXX))
allocate(cornerNodeTwoXX(numCornerPeriodicPairsXX))
allocate(cornerNodeThreeXX(numCornerPeriodicPairsXX))
allocate(cornerNodeFourXX(numCornerPeriodicPairsXX))

CALL nodePairingZZit%begin(nodePairingZZhash)

counter = 0
do keyIndex = 1, nodePairingZZhash%key_count()
  CALL nodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)

  srcZZ = nodePairingZZhashKey%ints(1)
  dstZZ = nodePairingZZvalue

  nodePairingXXhashInverseKey%ints(1) = dstZZ
  CALL nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, srcXX, success)

  !if ((success).and.(dstZZ.ne.dstTriple).and.((abs(nodeCoord(2,dstZZ)-boxLow(2))<tol).or.(abs(nodeCoord(2,dstZZ)-boxHigh(2))<tol))) then
  if (success.and.(dstZZ.ne.dstTriple)) then
    counter = counter + 1

    nodePairingXXhashInverseKey%ints(1) = srcZZ
    CALL nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, srcAux)

    cornerNodeOneXX(counter) = srcZZ
    cornerNodeTwoXX(counter) = srcAux

    cornerNodeThreeXX(counter) = dstZZ
    cornerNodeFourXX(counter)  = srcXX

    !write(456,*) "srcZZ: ", srcZZ, " -> dstZZ: ", dstZZ,    " [", nodeCoord(1,srcZZ), ", ", nodeCoord(2,srcZZ), ", ", nodeCoord(3,srcZZ), "] -> [", nodeCoord(1,dstZZ),   ", ", nodeCoord(2,dstZZ),   ", ", nodeCoord(3,dstZZ),   "]"
    !write(456,*) "dstXX: ", dstZZ,   " -> srcXX: ", srcXX , " [", nodeCoord(1,dstZZ),   ", ", nodeCoord(2,dstZZ),   ", ", nodeCoord(3,dstZZ),   "] -> [", nodeCoord(1,srcXX), ", ", nodeCoord(2,srcXX), ", ", nodeCoord(3,srcXX), "]"
    !write(456,*) "-----------------------"
  endif
enddo

deallocate(nodePairingYYhashInverseKey%ints)
deallocate(nodePairingXXhashInverseKey%ints)

return
!------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicCorners
