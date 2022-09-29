!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module geometry_mod
!----------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: nanopRadiusEff, numNanopFaces, wallDistance
use constants_mod,   only: pi
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
integer                              :: numDimensions, numNodes
integer                              :: numNodesLocalTypeFace, numNodesLocalTypeDomain
integer                              :: numElementsTypeFace, numElementsTypeDomain
integer                              :: numBulkNodePairs, numTotalNodePairs
integer                              :: numEdgePeriodicPairs
integer, allocatable, dimension(:)   :: nodePairId, numElementsOfNode
integer, allocatable, dimension(:)   :: nodeBelongsToFaceId
integer, allocatable, dimension(:)   :: edgeNodeOne, edgeNodeTwo, edgeNodeThree, edgeNodeFour
integer, allocatable, dimension(:,:) :: globalNodeIdTypeDomain, elementOfNode

logical, dimension(3,2)            :: isDirichletFace
logical, allocatable, dimension(:) :: nodeBelongsToDirichletFace
logical, allocatable, dimension(:) :: isDestPeriodicNodeXX, isDestPeriodicNodeYY, isDestPeriodicNodeZZ

type(fhash_type__ints_double) :: nodePairingXXhash, nodePairingYYhash, nodePairingZZhash
type(fhash_type__ints_double) :: nodePairingXXhashInverse, nodePairingYYhashInverse, nodePairingZZhashInverse

real(8), allocatable, dimension(:,:) :: nodeCoord
real(8), dimension(3)                :: boxLow, boxHigh, boxLength
real(8), parameter                   :: dx = 1.0d-1
real(8), parameter                   :: dy = 1.0d-1
real(8), parameter                   :: dz = 1.0d-1
!----------------------------------------------------------------------------------------------------------------------------!
  contains
    real(8) function interfaceArea()
      integer :: ii, mm, nn

      interfaceArea = 0.0d0

      do mm = 1, 3
        do nn = 1, 2
          if (isDirichletFace(mm,nn)) then
            interfaceArea = interfaceArea + DABS((boxHigh(1)-boxLow(1))*(boxHigh(2)-boxLow(2)))
          endif
        enddo
      enddo

      do ii = 1, numNanopFaces
        interfaceArea = interfaceArea + 4.0d0*pi*(nanopRadiusEff(ii)-wallDistance)**2.0d0
      enddo

      return
    end function interfaceArea
!----------------------------------------------------------------------------------------------------------------------------!
end module geometry_mod
