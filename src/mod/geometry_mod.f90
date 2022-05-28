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
integer                              :: nel, ndm, numnp, numel, nen_type_face, numel_type_face
integer                              :: numBulkNodePairs, numTotalNodePairs
integer, allocatable, dimension(:)   :: node_pair_id, num_of_elems_of_node
integer, allocatable, dimension(:)   :: nodeBelongsToFaceId
integer, allocatable, dimension(:,:) :: global_node_id_type_domain, el_node

logical, dimension(3,2)            :: isDirichletFace
logical, allocatable, dimension(:) :: nodeBelongsToDirichletFace

type(fhash_type__ints_double)          :: node_pairing_xx_hash, node_pairing_yy_hash, node_pairing_zz_hash, node_pairing_all_hash
type(fhash_type_iterator__ints_double) :: node_pairing_xx_it, node_pairing_yy_it, node_pairing_zz_it, node_pairing_all_it
type(ints_type)                        :: node_pairing_xx_key, node_pairing_yy_key, node_pairing_zz_key, node_pairing_all_key
integer                                :: node_pairing_xx_value, node_pairing_yy_value, node_pairing_zz_value, node_pairing_all_value

real(8), allocatable, dimension(:,:) :: xc
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
