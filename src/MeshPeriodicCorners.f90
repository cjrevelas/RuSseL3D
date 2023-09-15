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
                        destTriple, numCornerPeriodicPairsXX, &
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

integer :: sourceXX
integer :: sourceYY
integer :: sourceZZ, destZZ
integer :: sourceAux = 0
integer :: counter, keyIndex
!------------------------------------------------------------------------------------------------------!
allocate(nodePairingYYhashInverseKey%ints(1))
allocate(nodePairingXXhashInverseKey%ints(1))

CALL nodePairingZZit%begin(nodePairingZZhash)

counter = 0
do keyIndex = 1, nodePairingZZhash%key_count()
  CALL NodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)

  sourceZZ = nodePairingZZhashKey%ints(1)
  destZZ   = nodePairingZZvalue

  nodePairingYYhashInverseKey%ints(1) = destZZ
  CALL nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, sourceYY, success)

  !if ((success).and.((abs(nodeCoord(1,destZZ)-boxLow(1))<tol).or.(abs(nodeCoord(1,destZZ)-boxHigh(1))<tol))) counter = counter + 1
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

  sourceZZ = nodePairingZZhashKey%ints(1)
  destZZ   = nodePairingZZvalue

  nodePairingYYhashInverseKey%ints(1) = destZZ ! which is equal to destZZ
  CALL nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, sourceYY, success)

  !if ((success).and.((abs(nodeCoord(1,destZZ)-boxLow(1))<tol).or.(abs(nodeCoord(1,destZZ)-boxHigh(1))<tol))) then
  if (success) then
    counter = counter + 1

    nodePairingYYhashInverseKey%ints(1) = sourceZZ
    CALL nodePairingYYhashInverse%get(nodePairingYYhashInverseKey, sourceAux)

    cornerNodeOneYY(counter) = sourceZZ
    cornerNodeTwoYY(counter) = sourceAux

    cornerNodeThreeYY(counter) = destZZ
    cornerNodeFourYY(counter)  = sourceYY

    !write(123,*) "srcZZ: ", sourceZZ, " -> dstZZ: ", destZZ,    " [", nodeCoord(1,sourceZZ), ", ", nodeCoord(2,sourceZZ), ", ", nodeCoord(3,sourceZZ), "] -> [", nodeCoord(1,destZZ),   ", ", nodeCoord(2,destZZ),   ", ", nodeCoord(3,destZZ),   "]"
    !write(123,*) "dstYY: ", destZZ,   " -> srcYY: ", sourceYY , " [", nodeCoord(1,destZZ),   ", ", nodeCoord(2,destZZ),   ", ", nodeCoord(3,destZZ),   "] -> [", nodeCoord(1,sourceYY), ", ", nodeCoord(2,sourceYY), ", ", nodeCoord(3,sourceYY), "]"
    !write(123,*) "-----------------------"

    nodePairingXXhashInverseKey%ints(1) = destZZ
    CALL nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceXX, successAux)
    if (successAux) then
      destTriple = destZZ
      !write(123,*) "FOUND dstTriple: ", destZZ
    endif
  endif
enddo

CALL nodePairingZZit%begin(nodePairingZZhash)

counter = 0
do keyIndex = 1, nodePairingZZhash%key_count()
  CALL nodePairingZZit%next(nodePairingZZhashKey, nodePairingZZvalue)
  sourceZZ = nodePairingZZhashKey%ints(1)
  destZZ   = nodePairingZZvalue

  nodePairingXXhashInverseKey%ints(1) = destZZ
  CALL nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceXX, success)

  !if ((success).and.((abs(nodeCoord(2,destZZ)-boxLow(2))<tol).or.(abs(nodeCoord(2,destZZ)-boxHigh(2))<tol))) counter = counter + 1
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

  sourceZZ = nodePairingZZhashKey%ints(1)
  destZZ   = nodePairingZZvalue

  nodePairingXXhashInverseKey%ints(1) = destZZ
  CALL nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceXX, success)

  !if ((success).and.(destZZ.ne.destTriple).and.((abs(nodeCoord(2,destZZ)-boxLow(2))<tol).or.(abs(nodeCoord(2,destZZ)-boxHigh(2))<tol))) then
  if (success.and.(destZZ.ne.destTriple)) then
    counter = counter + 1

    nodePairingXXhashInverseKey%ints(1) = sourceZZ
    CALL nodePairingXXhashInverse%get(nodePairingXXhashInverseKey, sourceAux)

    cornerNodeOneXX(counter) = sourceZZ
    cornerNodeTwoXX(counter) = sourceAux

    cornerNodeThreeXX(counter) = destZZ
    cornerNodeFourXX(counter)  = sourceXX

    !write(456,*) "srcZZ: ", sourceZZ, " -> dstZZ: ", destZZ,    " [", nodeCoord(1,sourceZZ), ", ", nodeCoord(2,sourceZZ), ", ", nodeCoord(3,sourceZZ), "] -> [", nodeCoord(1,destZZ),   ", ", nodeCoord(2,destZZ),   ", ", nodeCoord(3,destZZ),   "]"
    !write(456,*) "dstXX: ", destZZ,   " -> srcXX: ", sourceXX , " [", nodeCoord(1,destZZ),   ", ", nodeCoord(2,destZZ),   ", ", nodeCoord(3,destZZ),   "] -> [", nodeCoord(1,sourceXX), ", ", nodeCoord(2,sourceXX), ", ", nodeCoord(3,sourceXX), "]"
    !write(456,*) "-----------------------"
  endif
enddo

deallocate(nodePairingYYhashInverseKey%ints)
deallocate(nodePairingXXhashInverseKey%ints)

return
!------------------------------------------------------------------------------------------------------!
end subroutine MeshPeriodicCorners
