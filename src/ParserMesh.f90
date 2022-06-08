!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ParserMesh()
!----------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use error_handing_mod
use write_helper_mod, only: adjl
use parser_vars_mod,  only: iow, periodicAxisId, domainIsPeriodic
use geometry_mod,     only: numNodesLocalTypeDomain, numElementsOfNode, boxLow, boxHigh, boxLength, &
                            nodeCoord, numNodes, numElementsTypeDomain, globalNodeIdTypeDomain,     &
                            numDimensions, nodePairId, numBulkNodePairs, numTotalNodePairs,         &
                            nodePairingXXhash, nodePairingYYhash, nodePairingZZhash,                &
                            numNodesLocalTypeFace, numElementsTypeFace
use kcw_mod,          only: F_m
use iofiles_mod,      only: IO_meshFile, IO_dirichletFaces, IO_nodepairs, IO_nodeConnectivity, &
                            IO_nodecoordinates, IO_meshProfile
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
character(len=200) :: line, aux_line

integer                              :: ii, jj, kk, pp
integer                              :: starting_pair, ending_pair
integer                              :: reason, id_of_entity_it_belongs, max_num_of_elems_per_node
integer                              :: numElementsTypeVertex, numElementsTypeEdge
integer                              :: numNodesLocalTypeVertex, numNodesLocalTypeEdge
integer, allocatable, dimension(:)   :: temp3
integer, allocatable, dimension(:,:) :: globalNodeIdTypeVertex, globalNodeIdTypeEdge, globalNodeIdTypeFace

real(8) :: box_volume = 0.0d0

logical :: success

type(fhash_type__ints_double) :: elemcon

type(fhash_type__ints_double)          :: vertex_entity_hash, edge_entity_hash, face_entity_hash, domain_entity_hash
type(fhash_type_iterator__ints_double) :: vertex_entity_it, edge_entity_it, face_entity_it, domain_entity_it
type(ints_type)                        :: vertex_entity_key, edge_entity_key, face_entity_key, domain_entity_key
integer                                :: vertex_entity_value, edge_entity_value, face_entity_value, domain_entity_value

type(fhash_type__ints_double)          :: xface1_hash, xface2_hash, yface1_hash, yface2_hash, zface1_hash, zface2_hash
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

      do ii = 1, numNodes
        read(12,*) (nodeCoord(jj,ii), jj = 1, numDimensions)

        do jj = 1, numDimensions
          boxHigh(jj) = max(nodeCoord(jj,ii), boxHigh(jj))
          boxLow(jj)  = min(nodeCoord(jj,ii), boxLow(jj))
        enddo
      enddo

      write(iow,*)
      write(6,*)

      write(iow,'(A6,A13,A17,A18)') "dim", "box_length", "box_min", "box_max"
      write(6  ,'(A6,A13,A17,A18)') "dim", "box_length", "box_min", "box_max"
      do jj = 1, numDimensions
        boxLength(jj) = boxHigh(jj) - boxLow(jj)
        write(iow,'(I5,2X,3(E16.9,2X))') jj, boxLength(jj), boxLow(jj), boxHigh(jj)
        write(6  ,'(I5,2X,3(E16.9,2X))') jj, boxLength(jj), boxLow(jj), boxHigh(jj)
      enddo

      write(6,*)
      write(iow,*)

      box_volume = boxLength(1) * boxLength(2) * boxLength(3)
      write(iow,'(3X,A40,E16.9,A13)')adjl("Box volume:",40), box_volume, " [Angstrom^3]"
      write(6  ,'(3X,A40,E16.9,A13)')adjl("Box volume:",40), box_volume, " [Angstrom^3]"
    elseif (INDEX(line,"3 vtx # type name")>0) then
      read(12,*)
      read(12,*)
      read(12,*) numNodesLocalTypeVertex
      read(12,*) numElementsTypeVertex
      allocate(globalNodeIdTypeVertex(numNodesLocalTypeVertex,numElementsTypeVertex))
      read(12,*)

      do kk = 1, numElementsTypeVertex
        read(12,*) (globalNodeIdTypeVertex(pp,kk), pp=1, numNodesLocalTypeVertex)
      enddo
      read(12,*)
      read(12,'(A100)',IOSTAT=reason) aux_line
      if (INDEX(aux_line," # number of parameter values per element")>0) then
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
      endif
      read(12,*)

      allocate(vertex_entity_key%ints(1))
      call vertex_entity_hash%reserve(numElementsTypeVertex)
      do kk = 1, numElementsTypeVertex
        vertex_entity_key%ints(1) = kk
        read(12,*) id_of_entity_it_belongs
        call vertex_entity_hash%get(vertex_entity_key, vertex_entity_value, success)
        if (.NOT.success) call vertex_entity_hash%set(vertex_entity_key, id_of_entity_it_belongs)
      enddo
    elseif ((INDEX(line,"3 edg # type name")>0).OR.(INDEX(line,"4 edg2 # type name")>0)) then
      read(12,*)
      read(12,*)
      read(12,*) numNodesLocalTypeEdge
      read(12,*) numElementsTypeEdge
      allocate(globalNodeIdTypeEdge(numNodesLocalTypeEdge,numElementsTypeEdge))
      read(12,*)

      do kk = 1, numElementsTypeEdge
        read(12,*) (globalNodeIdTypeEdge(pp,kk), pp=1,numNodesLocalTypeEdge)
      enddo
      read(12,*)
      read(12,'(A100)',IOSTAT=reason) aux_line
      if (INDEX(aux_line," # number of parameter values per element")>0) then
        read(12,*)
        read(12,*)
        do ii = 1, numElementsTypeEdge
          read(12,*)
        enddo
        read(12,*)
        read(12,*)
      endif
      read(12,*)

      allocate(edge_entity_key%ints(1))
      call edge_entity_hash%reserve(numElementsTypeEdge)
      do kk = 1, numElementsTypeEdge
        edge_entity_key%ints(1) = kk
        read(12,*) id_of_entity_it_belongs
        call edge_entity_hash%get(edge_entity_key, edge_entity_value, success)
        if (.NOT.success) call edge_entity_hash%set(edge_entity_key, id_of_entity_it_belongs)
      enddo
    elseif ((INDEX(line,"3 tri # type name")>0).OR.(INDEX(line,"4 tri2 # type name")>0)) then
      read(12,*)
      read(12,*)
      read(12,*) numNodesLocalTypeFace
      read(12,*) numElementsTypeFace
      read(12,*)

      allocate(globalNodeIdTypeFace(numNodesLocalTypeFace,numElementsTypeFace))
      do kk = 1, numElementsTypeFace
        read(12,*) (globalNodeIdTypeFace(pp,kk), pp=1,numNodesLocalTypeFace)
      enddo
      globalNodeIdTypeFace = globalNodeIdTypeFace + 1

      read(12,*)
      read(12,'(A100)',IOSTAT=reason) aux_line
      if ((INDEX(aux_line,"3 # number of parameter values per element")>0).OR.&
          (INDEX(aux_line,"6 # number of parameter values per element")>0)) then
        read(12,*)
        read(12,*)
        do ii = 1, numElementsTypeFace
          read(12,*)
        enddo
        read(12,*)
        read(12,*)
      endif
      read(12,*)

      allocate(face_entity_key%ints(1))
      call face_entity_hash%reserve(numElementsTypeFace)
      do kk = 1, numElementsTypeFace
        face_entity_key%ints(1) = kk
        read(12,*) id_of_entity_it_belongs
        call face_entity_hash%get(face_entity_key, face_entity_value, success)
        if (.NOT.success) call face_entity_hash%set(face_entity_key, id_of_entity_it_belongs+1)
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
      do ii = 1, numElementsTypeDomain
        read(12,*) (globalNodeIdTypeDomain(jj,ii), jj = 1, numNodesLocalTypeDomain)
      enddo
      globalNodeIdTypeDomain = globalNodeIdTypeDomain + 1

      if (numNodesLocalTypeDomain>4) then
        allocate(temp3(numElementsTypeDomain))
        do ii = 1, numElementsTypeDomain
          temp3(ii) = globalNodeIdTypeDomain(7,ii)
          globalNodeIdTypeDomain(7,ii) = globalNodeIdTypeDomain(6,ii)
          globalNodeIdTypeDomain(6,ii) = temp3(ii)
        enddo
      endif

      read(12,*)
      read(12,'(A100)',IOSTAT=reason) aux_line
      if ((INDEX(aux_line,"4 # number of parameter values per element")>0).OR.(INDEX(aux_line,"10 # number of parameter values per element")>0)) then
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
      endif
      read(12,*)

      allocate(domain_entity_key%ints(1))
      call domain_entity_hash%reserve(numElementsTypeDomain)
      do kk = 1, numElementsTypeDomain
        domain_entity_key%ints(1) = kk
        read(12,*) id_of_entity_it_belongs
        call domain_entity_hash%get(domain_entity_key, domain_entity_value, success)
        if (.NOT.success) call domain_entity_hash%set(domain_entity_key, id_of_entity_it_belongs)
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

call vertex_entity_it%begin(vertex_entity_hash)
call edge_entity_it%begin(edge_entity_hash)
call domain_entity_it%begin(domain_entity_hash)

do kk = 1, vertex_entity_hash%key_count()
  call vertex_entity_it%next(vertex_entity_key, vertex_entity_value)
  write(11,*) vertex_entity_value, vertex_entity_key%ints(1)
enddo

do kk = 1, edge_entity_hash%key_count()
  call edge_entity_it%next(edge_entity_key, edge_entity_value)
  write(22,*) edge_entity_value, edge_entity_key%ints(1)
enddo

do kk = 1, domain_entity_hash%key_count()
  call domain_entity_it%next(domain_entity_key, domain_entity_value)
  write(44,*) domain_entity_value, domain_entity_key%ints(1)
enddo
#endif

if (domainIsPeriodic) then
  if (periodicAxisId(1)) call MeshFaceEntities('x', face_entity_hash, xface1_hash, xface2_hash)
  if (periodicAxisId(2)) call MeshFaceEntities('y', face_entity_hash, yface1_hash, yface2_hash)
  if (periodicAxisId(3)) call MeshFaceEntities('z', face_entity_hash, zface1_hash, zface2_hash)

#ifdef DEBUG_OUTPUTS
  call face_entity_it%begin(face_entity_hash)
  do kk = 1, face_entity_hash%key_count()
    call face_entity_it%next(face_entity_key, face_entity_value)
    write(33,*) face_entity_value-1, face_entity_key%ints(1)
  enddo
  close(11)
  close(22)
  close(33)
  close(44)

  if (periodicAxisId(1)) call MeshPeriodicFaces('x', xface1_hash, xface2_hash)
  if (periodicAxisId(2)) call MeshPeriodicFaces('y', yface1_hash, yface2_hash)
  if (periodicAxisId(3)) call MeshPeriodicFaces('z', zface1_hash, zface2_hash)
#endif

  if (periodicAxisId(1)) call MeshPeriodicNodePairs(globalNodeIdTypeFace, 'x', xface1_hash, xface2_hash, nodePairingXXhash)
  if (periodicAxisId(2)) call MeshPeriodicNodePairs(globalNodeIdTypeFace, 'y', yface1_hash, yface2_hash, nodePairingYYhash)
  if (periodicAxisId(3)) call MeshPeriodicNodePairs(globalNodeIdTypeFace, 'z', zface1_hash, zface2_hash, nodePairingZZhash)
endif

allocate(numElementsOfNode(numNodes))
call MeshElementsOfNode(max_num_of_elems_per_node, numElementsOfNode)

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

call MeshBulkNodePairs(nodePairingXXhash, nodePairingYYhash, nodePairingZZhash, elemcon)

! xx pairs
starting_pair = numBulkNodePairs + 1
ending_pair   = numBulkNodePairs + nodePairingXXhash%key_count()

if (periodicAxisId(1)) call MeshAppendPeriodicPairs(elemcon, starting_pair, ending_pair, nodePairingXXhash)

do ii = 1, numBulkNodePairs + nodePairingXXhash%key_count()
  F_m%is_zero(ii) = (nodePairId(ii)/=ii)
enddo

! yy pairs
starting_pair = numBulkNodePairs + nodePairingXXhash%key_count() + 1
ending_pair   = numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count()

if (periodicAxisId(2)) call MeshAppendPeriodicPairs(elemcon, starting_pair, ending_pair, nodePairingYYhash)

do ii = numBulkNodePairs + nodePairingXXhash%key_count() + 1, &
        numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count()
  F_m%is_zero(ii) = (nodePairId(ii)/=ii)
enddo

! zz pairs
starting_pair = numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count() + 1
ending_pair   = numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count() + nodePairingZZhash%key_count()

if (periodicAxisId(3)) call MeshAppendPeriodicPairs(elemcon, starting_pair, ending_pair, nodePairingZZhash)

do ii = numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count() + 1, &
        numBulkNodePairs + nodePairingXXhash%key_count() + nodePairingYYhash%key_count() + nodePairingZZhash%key_count()
  F_m%is_zero(ii) = (nodePairId(ii)/=ii)
enddo

call elemcon%clear()

call MeshDirichletFaces(numElementsTypeFace, numNodesLocalTypeFace, globalNodeIdTypeFace, face_entity_hash)

#ifdef DEBUG_OUTPUTS
open(unit=77, file = IO_nodePairs)
do ii = 1, numTotalNodePairs
  write(77,'(4(2X,I9),2X,L9)') ii, F_m%row(ii), F_m%col(ii), nodePairId(ii), F_m%is_zero(ii)
enddo
close(77)

open (unit=77, file = IO_nodeConnectivity)
do ii = 1, numElementsTypeDomain
  write(77,'(11(2X,I9))') (globalNodeIdTypeDomain(jj,ii), jj = 1, numNodesLocalTypeDomain)
enddo
close(77)

open(77, file = IO_nodeCoordinates)
write(77,'(3(2X,A16))') 'x', 'y', 'z'
do ii = 1, numNodes
  write(77,'(3(2X,F16.9))') (nodeCoord(jj,ii), jj = 1, numDimensions)
enddo
close(77)

call MeshProfile()
#endif

#ifdef DEBUG_OUTPUTS
deallocate(vertex_entity_key%ints, edge_entity_key%ints, domain_entity_key%ints)
call vertex_entity_hash%clear()
call edge_entity_hash%clear()
call domain_entity_hash%clear()
#endif

! Deallocate memory
if (numNodesLocalTypeDomain > 4) deallocate(temp3)
deallocate(globalNodeIdTypeVertex, globalNodeIdTypeEdge, globalNodeIdTypeFace)

return
!----------------------------------------------------------------------------------------------------------------------------!
end subroutine ParserMesh
