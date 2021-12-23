!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module geometry_mod

  use, intrinsic :: iso_fortran_env
  use fhash_module__ints_double
  use ints_module
  use parser_vars_mod, only: radius_np_eff, num_of_nanoparticle_faces, wall_distance
  use constants_mod,   only: pi

  implicit none

    integer                              :: nel, ndm, numnp, numel
    integer                              :: num_of_bulk_pairs, total_num_of_node_pairs
    integer                              :: num_dest_xx_neighbors, num_dest_yy_neighbors, num_dest_zz_neighbors
    integer, allocatable, dimension(:)   :: node_pair_id, num_of_elems_of_node
    integer, allocatable, dimension(:,:) :: global_node_id_type_domain, el_node

    logical, dimension(3,2)            :: is_dirichlet_face
    logical, allocatable, dimension(:) :: node_belongs_to_dirichlet_face

    type(fhash_type__ints_double)          :: node_pairing_xx_hash, node_pairing_yy_hash, node_pairing_zz_hash, node_pairing_all_hash
    type(fhash_type_iterator__ints_double) :: node_pairing_xx_it, node_pairing_yy_it, node_pairing_zz_it, node_pairing_all_it
    type(ints_type)                        :: node_pairing_xx_key, node_pairing_yy_key, node_pairing_zz_key, node_pairing_all_key
    integer                                :: node_pairing_xx_value, node_pairing_yy_value, node_pairing_zz_value, node_pairing_all_value

    real(8), allocatable, dimension(:,:) :: xc
    real(8), dimension(3)                :: box_lo, box_hi, box_len
    real(8), parameter                   :: dx = 1.0d-1
    real(8), parameter                   :: dy = 1.0d-1
    real(8), parameter                   :: dz = 1.0d-1

  contains

    real(8) function interf_area()
        integer :: ii, mm, nn

        interf_area = 0.0d0

        do mm = 1, 3
            do nn = 1, 2
                if (is_dirichlet_face(mm,nn)) then
                    interf_area = interf_area + DABS((box_hi(1)-box_lo(1))*(box_hi(2)-box_lo(2)))
                endif
            enddo
        enddo

        do ii = 1, num_of_nanoparticle_faces
            interf_area = interf_area + 4.0d0*pi*(radius_np_eff(ii)-wall_distance)**2
        enddo

        return
    end function interf_area
end module geometry_mod
