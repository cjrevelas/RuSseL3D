!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module geometry_mod
!----------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use parser_vars_mod, only: radius_np_eff, numNanoparticleFaces, wall_distance
use constants_mod,   only: pi
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
integer                              :: numDimensions, numNodes
integer                              :: numNodesLocalTypeFace, numNodesLocalTypeDomain
integer                              :: numElementsTypeFace, numElementsTypeDomain
integer                              :: numBulkNodePairs, numTotalNodePairs
integer, allocatable, dimension(:)   :: nodePairId, numElementsOfNode
integer, allocatable, dimension(:)   :: nodeBelongsToFaceId
integer, allocatable, dimension(:,:) :: globalNodeIdTypeDomain, elementOfNode

logical, dimension(3,2)            :: isDirichletFace
logical, allocatable, dimension(:) :: nodeBelongsToDirichletFace

type(fhash_type__ints_double) :: nodePairingXXhash, nodePairingYYhash, nodePairingZZhash

real(8), allocatable, dimension(:,:) :: nodeCoord
real(8), dimension(3)                :: boxLow, boxHigh, boxLength
real(8), parameter                   :: dx = 1.0d-1
real(8), parameter                   :: dy = 1.0d-1
real(8), parameter                   :: dz = 1.0d-1
!----------------------------------------------------------------------------------------------------------------------------!
  contains

    real(8) function interf_area()
      integer :: ii, mm, nn

      interf_area = 0.0d0

      do mm = 1, 3
        do nn = 1, 2
          if (isDirichletFace(mm,nn)) then
            interf_area = interf_area + DABS((boxHigh(1)-boxLow(1))*(boxHigh(2)-boxLow(2)))
          endif
        enddo
      enddo

      do ii = 1, numNanoparticleFaces
        interf_area = interf_area + 4.0d0*pi*(radius_np_eff(ii)-wall_distance)**2.0d0
      enddo

      return
    end function interf_area
!----------------------------------------------------------------------------------------------------------------------------!
end module geometry_mod
