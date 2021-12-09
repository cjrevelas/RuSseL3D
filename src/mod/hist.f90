module hist
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer                                :: nbin
integer, allocatable, dimension(:,:,:) :: planar_cell_of_np
integer, allocatable, dimension(:,:)   :: sph_cell_of_np

real(8)                                :: lbin
real(8), allocatable, dimension(:,:,:) :: dist_from_face, cell_vol_planar
real(8), allocatable, dimension(:,:)   :: dist_from_np, cell_vol_sph
!--------------------------------------------------------------------!
end module hist
