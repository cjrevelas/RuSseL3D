!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module hist_mod
    integer                                :: nbin
    integer, allocatable, dimension(:,:,:) :: planar_cell_of_np
    integer, allocatable, dimension(:,:)   :: sph_cell_of_np

    real(8)                                :: lbin
    real(8), allocatable, dimension(:,:,:) :: dist_from_face, cell_vol_planar
    real(8), allocatable, dimension(:,:)   :: dist_from_np, cell_vol_sph
end module hist_mod
