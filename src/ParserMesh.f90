!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ParserMesh(elemcon)
!----------------------------------------------------------------------------------------------------------------------------!
use fhash_module__ints_double
use ints_module
use error_handing_mod
use write_helper_mod, only: adjl
use parser_vars_mod,  only: iow, periodicAxisId, domainIsPeriodic
use geometry_mod,     only: numNodesLocalTypeDomain, numElementsOfNode, boxLow, boxHigh, boxLength,       &
                            nodeCoord, numNodes, numElementsTypeDomain, globalNodeIdTypeDomain,           &
                            numDimensions, nodePairId, numBulkNodePairs, numTotalNodePairs,               &
                            nodePairingXXhash, nodePairingYYhash, nodePairingZZhash,                      &
                            nodePairingXXhashInverse, nodePairingYYhashInverse, nodePairingZZhashInverse, &
                            numNodesLocalTypeFace, numElementsTypeFace, isDestPeriodicNodeXX, &
                            isDestPeriodicNodeYY, isDestPeriodicNodeZZ
use kcw_mod,          only: F_m
use iofiles_mod,      only: IO_meshFile, IO_dirichletFaces, IO_nodepairs, IO_nodeConnectivity, &
                            IO_nodecoordinates, IO_meshProfile
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
character(len=200) :: line, auxLine

integer                              :: element, localNodeId, globalNodeId, pairId, axis, keyIndex
integer                              :: startingPair, endingPair
integer                              :: reason, idOfEntityItBelongs, maxNumOfElementsPerNode
integer                              :: numElementsTypeVertex, numElementsTypeEdge
integer                              :: numNodesLocalTypeVertex, numNodesLocalTypeEdge
integer, allocatable, dimension(:)   :: temp3
integer, allocatable, dimension(:,:) :: globalNodeIdTypeVertex, globalNodeIdTypeEdge, globalNodeIdTypeFace

real(8) :: boxVolume = 0.0d0

logical :: success

type(fhash_type__ints_double), intent(inout) :: elemcon

type(fhash_type__ints_double)          :: vertexEntityHash, edgeEntityHash, faceEntityHash, domainEntityHash
type(fhash_type_iterator__ints_double) :: vertexEntityIt, edgeEntityIt, faceEntityIt, domainEntityIt
type(ints_type)                        :: vertexEntityKey, edgeEntityKey, faceEntityKey, domainEntityKey
integer                                :: vertexEntityValue, edgeEntityValue, faceEntityValue, domainEntityValue

type(fhash_type__ints_double)          :: xfaceOneHash, xfaceTwoHash, yfaceOneHash, yfaceTwoHash, zfaceOneHash, zfaceTwoHash
!----------------------------------------------------------------------------------------------------------------------------!
open(unit=12, file=IO_meshFile)

do
  read(12,'(A100)',IOSTAT=reason) line

  if (reason>0) then
    write(*,*) "Something went wrong!"
  elseif (reason<0) then
    exit
  else
    if (INDEX(line,"# sdim") > 0) then
      read(line,*) numDimensions
    elseif (INDEX(line,"# number of mesh points")>0) then
      read(line,*) numNodes
      allocate(nodeCoord(numDimensions,numNodes))
    elseif (INDEX(line,"# Mesh point coordinates")>0) then
      boxLow    = 0.0d0
      boxHigh   = 0.0d0
      boxLength = 0.0d0

      do globalNodeId = 1, numNodes
        read(12,*) (nodeCoord(axis,globalNodeId), axis = 1, numDimensions)

        do axis = 1, numDimensions
          boxHigh(axis) = MAX(nodeCoord(axis,globalNodeId), boxHigh(axis))
          boxLow(axis)  = MIN(nodeCoord(axis,globalNodeId), boxLow(axis))
        enddo
      enddo

      write(iow,*)
      write(6,*)

      write(iow,'(A6,A13,A17,A18)') "dim", "box_length", "box_min", "box_max"
      write(6  ,'(A6,A13,A17,A18)') "dim", "box_length", "box_min", "box_max"
      do axis = 1, numDimensions
        boxLength(axis) = boxHigh(axis) - boxLow(axis)
        write(iow,'(I5,2X,3(E16.9,2X))') axis, boxLength(axis), boxLow(axis), boxHigh(axis)
        write(6  ,'(I5,2X,3(E16.9,2X))') axis, boxLength(axis), boxLow(axis), boxHigh(axis)
      enddo

      write(6,*)
      write(iow,*)

      boxVolume = boxLength(1) * boxLength(2) * boxLength(3)
      write(iow,'(3X,A40,E16.9,A13)')adjl("Box volume:",40), boxVolume, " [Angstrom^3]"
      write(6  ,'(3X,A40,E16.9,A13)')adjl("Box volume:",40), boxVolume, " [Angstrom^3]"
    elseif (INDEX(line,"3 vtx # type name")>0) then
      read(12,*)
      read(12,*)
      read(12,*) numNodesLocalTypeVertex
      read(12,*) numElementsTypeVertex
      allocate(globalNodeIdTypeVertex(numNodesLocalTypeVertex,numElementsTypeVertex))
      read(12,*)

      do element = 1, numElementsTypeVertex
        read(12,*) (globalNodeIdTypeVertex(localNodeId,element), localNodeId = 1, numNodesLocalTypeVertex)
      enddo
      read(12,*)
      read(12,'(A100)',IOSTAT=reason) auxLine
      if (INDEX(auxLine," # number of parameter values per element")>0) then
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
      endif
      read(12,*)

      allocate(vertexEntityKey%ints(1))
      call vertexEntityHash%reserve(numElementsTypeVertex)
      do element = 1, numElementsTypeVertex
        vertexEntityKey%ints(1) = element
        read(12,*) idOfEntityItBelongs
        call vertexEntityHash%get(vertexEntityKey, vertexEntityValue, success)
        if (.NOT.success) call vertexEntityHash%set(vertexEntityKey, idOfEntityItBelongs)
      enddo
    elseif ((INDEX(line,"3 edg # type name")>0).OR.(INDEX(line,"4 edg2 # type name")>0)) then
      read(12,*)
      read(12,*)
      read(12,*) numNodesLocalTypeEdge
      read(12,*) numElementsTypeEdge
      allocate(globalNodeIdTypeEdge(numNodesLocalTypeEdge,numElementsTypeEdge))
      read(12,*)

      do element = 1, numElementsTypeEdge
        read(12,*) (globalNodeIdTypeEdge(localNodeId,element), localNodeId = 1, numNodesLocalTypeEdge)
      enddo
      read(12,*)
      read(12,'(A100)',IOSTAT=reason) auxLine
      if (INDEX(auxLine," # number of parameter values per element")>0) then
        read(12,*)
        read(12,*)
        do element = 1, numElementsTypeEdge
          read(12,*)
        enddo
        read(12,*)
        read(12,*)
      endif
      read(12,*)

      allocate(edgeEntityKey%ints(1))
      call edgeEntityHash%reserve(numElementsTypeEdge)
      do element = 1, numElementsTypeEdge
        edgeEntityKey%ints(1) = element
        read(12,*) idOfEntityItBelongs
        call edgeEntityHash%get(edgeEntityKey, edgeEntityValue, success)
        if (.NOT.success) call edgeEntityHash%set(edgeEntityKey, idOfEntityItBelongs)
      enddo
    elseif ((INDEX(line,"3 tri # type name")>0).OR.(INDEX(line,"4 tri2 # type name")>0)) then
      read(12,*)
      read(12,*)
      read(12,*) numNodesLocalTypeFace
      read(12,*) numElementsTypeFace
      read(12,*)

      allocate(globalNodeIdTypeFace(numNodesLocalTypeFace,numElementsTypeFace))
      do element = 1, numElementsTypeFace
        read(12,*) (globalNodeIdTypeFace(localNodeId,element), localNodeId = 1, numNodesLocalTypeFace)
      enddo
      globalNodeIdTypeFace = globalNodeIdTypeFace + 1

      read(12,*)
      read(12,'(A100)',IOSTAT=reason) auxLine
      if ((INDEX(auxLine,"3 # number of parameter values per element")>0).OR.&
          (INDEX(auxLine,"6 # number of parameter values per element")>0)) then
        read(12,*)
        read(12,*)
        do element = 1, numElementsTypeFace
          read(12,*)
        enddo
        read(12,*)
        read(12,*)
      endif
      read(12,*)

      allocate(faceEntityKey%ints(1))
      call faceEntityHash%reserve(numElementsTypeFace)
      do element = 1, numElementsTypeFace
        faceEntityKey%ints(1) = element
        read(12,*) idOfEntityItBelongs
        call faceEntityHash%get(faceEntityKey, faceEntityValue, success)
        if (.NOT.success) call faceEntityHash%set(faceEntityKey, idOfEntityItBelongs+1)
      enddo
    elseif ((INDEX(line,"3 tet # type name")>0).OR.(INDEX(line,"4 tet2 # type name")>0)) then
      read(12,*)
      read(12,*)
      read(12,*) numNodesLocalTypeDomain
      read(12,*) numElementsTypeDomain
      read(12,*)

      numBulkNodePairs = numNodesLocalTypeDomain * numNodesLocalTypeDomain * numElementsTypeDomain

      write(iow,*)
      write(*,*)
      write(iow,'(3X,"Number of mesh points:                 ",I16)') numNodes
      write(iow,'(3X,"Number of elements:                    ",I16)') numElementsTypeDomain
      write(iow,'(3X,"Number of nodes per element:           ",I16)') numNodesLocalTypeDomain
      write(iow,'(3X,"Number of matrix indeces:              ",I16)') numBulkNodePairs
      write(6,'(3X,"Number of mesh points                    ",I16)') numNodes
      write(6,'(3X,"Number of elements:                      ",I16)') numElementsTypeDomain
      write(6,'(3X,"Number of nodes per element:             ",I16)') numNodesLocalTypeDomain
      write(6,'(3X,"Number of matrix indeces:                ",I16)') numBulkNodePairs

      allocate(globalNodeIdTypeDomain(numNodesLocalTypeDomain,numElementsTypeDomain))
      do element = 1, numElementsTypeDomain
        read(12,*) (globalNodeIdTypeDomain(localNodeId,element), localNodeId = 1, numNodesLocalTypeDomain)
      enddo
      globalNodeIdTypeDomain = globalNodeIdTypeDomain + 1

      if (numNodesLocalTypeDomain>4) then
        allocate(temp3(numElementsTypeDomain))
        do element = 1, numElementsTypeDomain
          temp3(element) = globalNodeIdTypeDomain(7,element)
          globalNodeIdTypeDomain(7,element) = globalNodeIdTypeDomain(6,element)
          globalNodeIdTypeDomain(6,element) = temp3(element)
        enddo
      endif

      read(12,*)
      read(12,'(A100)',IOSTAT=reason) auxLine
      if ((INDEX(auxLine,"4 # number of parameter values per element")>0).OR.(INDEX(auxLine,"10 # number of parameter values per element")>0)) then
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
      endif
      read(12,*)

      allocate(domainEntityKey%ints(1))
      call domainEntityHash%reserve(numElementsTypeDomain)
      do element = 1, numElementsTypeDomain
        domainEntityKey%ints(1) = element
        read(12,*) idOfEntityItBelongs
        call domainEntityHash%get(domainEntityKey, domainEntityValue, success)
        if (.NOT.success) call domainEntityHash%set(domainEntityKey, idOfEntityItBelongs)
      enddo
    endif
  endif
enddo
close(12)

#ifdef DEBUG_OUTPUTS
open(unit=11, file=IO_vertexElements)
open(unit=22, file=IO_edgeElements)
open(unit=33, file=IO_faceElements)
open(unit=44, file=IO_domainElements)

call vertexEntityIt%begin(vertexEntityHash)
call edgeEntityIt%begin(edgeEntityHash)
call domainEntityIt%begin(domainEntityHash)

do keyIndex = 1, vertexEntityHash%key_count()
  call vertexEntityIt%next(vertexEntityKey, vertexEntityValue)
  write(11,*) vertexEntityValue, vertexEntityKey%ints(1)
enddo

do keyIndex = 1, edgeEntityHash%key_count()
  call edgeEntityIt%next(edgeEntityKey, edgeEntityValue)
  write(22,*) edgeEntityValue, edgeEntityKey%ints(1)
enddo

do keyIndex = 1, domainEntityHash%key_count()
  call domainEntityIt%next(domainEntityKey, domainEntityValue)
  write(44,*) domainEntityValue, domainEntityKey%ints(1)
enddo
#endif

if (domainIsPeriodic) then
  if (periodicAxisId(1)) then
    call MeshFaceEntities(1, boxLow(1),  globalNodeIdTypeFace, faceEntityHash, xfaceOneHash)
    call MeshFaceEntities(1, boxHigh(1), globalNodeIdTypeFace, faceEntityHash, xfaceTwoHash)
  endif
  if (periodicAxisId(2)) then
    call MeshFaceEntities(2, boxLow(2),  globalNodeIdTypeFace, faceEntityHash, yfaceOneHash)
    call MeshFaceEntities(2, boxHigh(2), globalNodeIdTypeFace, faceEntityHash, yfaceTwoHash)
  endif
  if (periodicAxisId(3)) then
    call MeshFaceEntities(3, boxHigh(3),  globalNodeIdTypeFace, faceEntityHash, zfaceOneHash)
    call MeshFaceEntities(3, boxLow(3),   globalNodeIdTypeFace, faceEntityHash, zfaceTwoHash)
  endif

#ifdef DEBUG_OUTPUTS
  call faceEntityIt%begin(faceEntityHash)
  do keyIndex = 1, faceEntityHash%key_count()
    call faceEntityIt%next(faceEntityKey, faceEntityValue)
    write(33,*) faceEntityValue-1, faceEntityKey%ints(1)
  enddo
  close(11)
  close(22)
  close(33)
  close(44)

  if (periodicAxisId(1)) call MeshPeriodicFaces('x', xfaceOneHash, xfaceTwoHash)
  if (periodicAxisId(2)) call MeshPeriodicFaces('y', yfaceOneHash, yfaceTwoHash)
  if (periodicAxisId(3)) call MeshPeriodicFaces('z', zfaceOneHash, zfaceTwoHash)
#endif

  if (periodicAxisId(1)) then
    allocate(isDestPeriodicNodeXX(numNodes))
    call MeshPeriodicNodePairs(globalNodeIdTypeFace, 'x', xfaceOneHash, xfaceTwoHash, isDestPeriodicNodeXX, nodePairingXXhash, nodePairingXXhashInverse)
  endif
  if (periodicAxisId(2)) then
    allocate(isDestPeriodicNodeYY(numNodes))
    call MeshPeriodicNodePairs(globalNodeIdTypeFace, 'y', yfaceOneHash, yfaceTwoHash, isDestPeriodicNodeYY, nodePairingYYhash, nodePairingYYhashInverse)
  endif
  if (periodicAxisId(3)) then
    allocate(isDestPeriodicNodeZZ(numNodes))
    call MeshPeriodicNodePairs(globalNodeIdTypeFace, 'z', zfaceOneHash, zfaceTwoHash, isDestPeriodicNodeZZ, nodePairingZZhash, nodePairingZZhashInverse)
  endif
endif

allocate(numElementsOfNode(numNodes))
call MeshElementsOfNode(maxNumOfElementsPerNode, numElementsOfNode)

numTotalNodePairs = numBulkNodePairs

if (periodicAxisId(1)) numTotalNodePairs = numTotalNodePairs + nodePairingXXhash%key_count()
if (periodicAxisId(2)) numTotalNodePairs = numTotalNodePairs + nodePairingYYhash%key_count()
if (periodicAxisId(3)) numTotalNodePairs = numTotalNodePairs + nodePairingZZhash%key_count()

allocate(F_m%row(numTotalNodePairs))
allocate(F_m%col(numTotalNodePairs))
allocate(F_m%g(numTotalNodePairs))
allocate(F_m%rh(numTotalNodePairs))
allocate(F_m%c(numTotalNodePairs))
allocate(F_m%k(numTotalNodePairs))
allocate(F_m%w(numTotalNodePairs))
allocate(F_m%is_zero(numTotalNodePairs))

F_m%g       = 0.0d0
F_m%k       = 0.0d0
F_m%c       = 0.0d0
F_m%rh      = 0.0d0
F_m%row     = 0
F_m%col     = 0
F_m%is_zero = .True.

call MeshBulkNodePairs(elemcon)

if (domainIsPeriodic) then
  call MeshPeriodicDestNeighbors(isDestPeriodicNodeXX, isDestPeriodicNodeYY, isDestPeriodicNodeZZ, elemcon, nodePairId)
  call MeshPeriodicEdges(nodePairingXXhash, nodePairingYYhash, nodePairingYYhashInverse)
  call MeshPeriodicCorners()
endif

! xx pairs
startingPair = numBulkNodePairs + 1
endingPair   = numBulkNodePairs + nodePairingXXhash%key_count()

if (periodicAxisId(1)) call MeshAppendPeriodicPairs(elemcon, startingPair, endingPair, nodePairingXXhash)

do pairId = 1, numBulkNodePairs + nodePairingXXhash%key_count()
  F_m%is_zero(pairId) = (nodePairId(pairId) /= pairId)
enddo

! yy pairs
startingPair = numBulkNodePairs + nodePairingXXhash%key_count() + 1
endingPair   = numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count()

if (periodicAxisId(2)) call MeshAppendPeriodicPairs(elemcon, startingPair, endingPair, nodePairingYYhash)

do pairId = numBulkNodePairs + nodePairingXXhash%key_count() + 1, &
        numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count()
  F_m%is_zero(pairId) = (nodePairId(pairId) /= pairId)
enddo

! zz pairs
startingPair = numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count() + 1
endingPair   = numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count() + nodePairingZZhash%key_count()

if (periodicAxisId(3)) call MeshAppendPeriodicPairs(elemcon, startingPair, endingPair, nodePairingZZhash)

do pairId = numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count() + 1, &
        numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count() + nodePairingZZhash%key_count()
  F_m%is_zero(pairId) = (nodePairId(pairId) /= pairId)
enddo

call MeshDirichletFaces(numElementsTypeFace, numNodesLocalTypeFace, globalNodeIdTypeFace, faceEntityHash)

#ifdef DEBUG_OUTPUTS
open(unit=77, file = IO_nodePairs)
do pairId = 1, numTotalNodePairs
  write(77,'(4(2X,I9),2X,L9)') pairId, F_m%row(pairId), F_m%col(pairId), nodePairId(pairId), F_m%is_zero(pairId)
enddo
close(77)

open (unit=77, file = IO_nodeConnectivity)
do element = 1, numElementsTypeDomain
  write(77,'(11(2X,I9))') (globalNodeIdTypeDomain(localNodeId,element), localNodeId = 1, numNodesLocalTypeDomain)
enddo
close(77)

open(77, file = IO_nodeCoordinates)
write(77,'(3(2X,A16))') 'x', 'y', 'z'
do element = 1, numNodes
  write(77,'(3(2X,F16.9))') (nodeCoord(axis,element), axis = 1, numDimensions)
enddo
close(77)

call MeshProfile()
#endif

#ifdef DEBUG_OUTPUTS
deallocate(vertexEntityKey%ints, edgeEntityKey%ints, domainEntityKey%ints)
call vertexEntityHash%clear()
call edgeEntityHash%clear()
call domainEntityHash%clear()
#endif

! Deallocate memory
if (numNodesLocalTypeDomain > 4) deallocate(temp3)
deallocate(globalNodeIdTypeVertex, globalNodeIdTypeEdge, globalNodeIdTypeFace)
if (periodicAxisId(1)) deallocate(isDestPeriodicNodeXX)
if (periodicAxisId(2)) deallocate(isDestPeriodicNodeYY)
if (periodicAxisId(3)) deallocate(isDestPeriodicNodeZZ)

return
!----------------------------------------------------------------------------------------------------------------------------!
end subroutine ParserMesh
