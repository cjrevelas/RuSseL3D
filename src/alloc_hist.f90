subroutine alloc_hist(upd_lbin, volnp)
!-----------------------------------------------------------------------------------------------------------!
use hist
use geometry,    only: numnp, box_len, is_dir_face, box_lo, box_hi, xc
use parser_vars, only: n_nanopart_faces, wall_distance, center_np, radius_np_eff
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: kk, mm, nn, bin

real(8), intent(in)                   :: upd_lbin
real(8), intent(in), dimension(numnp) :: volnp
real(8)                               :: r_center_surf, r_centers, radius_np_actual, Rmax
!-----------------------------------------------------------------------------------------------------------!
lbin = upd_lbin
Rmax = SQRT(box_len(1)**2 + box_len(2)**2 + box_len(3)**2)
nbin = INT(Rmax / lbin) + 1

if (ALLOCATED(planar_cell_of_np)) deallocate(planar_cell_of_np)
if (ALLOCATED(sph_cell_of_np))    deallocate(sph_cell_of_np)
if (ALLOCATED(dist_from_face))    deallocate(dist_from_face)
if (ALLOCATED(dist_from_np))      deallocate(dist_from_np)
if (ALLOCATED(cell_vol_planar))   deallocate(cell_vol_planar)
if (ALLOCATED(cell_vol_sph))      deallocate(cell_vol_sph)

allocate(planar_cell_of_np(numnp,3,2))
allocate(dist_from_face(numnp,3,2))
allocate(cell_vol_planar(nbin,3,2))
allocate(sph_cell_of_np(n_nanopart_faces,numnp))
allocate(dist_from_np(n_nanopart_faces,numnp))
allocate(cell_vol_sph(n_nanopart_faces,nbin))

planar_cell_of_np = 0
sph_cell_of_np    = 0
dist_from_face    = 0.d0
dist_from_np      = 0.d0
cell_vol_planar   = 0.d0
cell_vol_sph      = 0.d0

! chains grafted to planar surfaces
do mm = 1, 3
    do nn = 1, 2
        if (is_dir_face(mm,nn)) then
            do kk = 1, numnp
                if (nn.eq.1) then
                    r_center_surf = xc(mm,kk) - box_lo(mm) + wall_distance
                elseif (nn.eq.2) then 
                    r_center_surf = box_hi(mm) - xc(mm,kk) + wall_distance
                endif
                bin                         = INT(r_center_surf/lbin)+1
                planar_cell_of_np(kk,mm,nn) = bin
                cell_vol_planar(bin,mm,nn)  = cell_vol_planar(bin,mm,nn) + volnp(kk)
                dist_from_face(kk,mm,nn)    = r_center_surf
            enddo
        endif
    enddo
enddo

! chains grafted to nanoparticles
do mm = 1, n_nanopart_faces
    do kk = 1, numnp
        r_centers             = DSQRT(  (xc(1,kk)-center_np(1,mm))**2 + (xc(2,kk)-center_np(2,mm))**2 &
&                             + (xc(3,kk)-center_np(3,mm))**2) 
        radius_np_actual      = radius_np_eff(mm) - wall_distance
        r_center_surf         = r_centers - radius_np_actual
        bin                   = INT(r_center_surf/lbin)+1
        sph_cell_of_np(mm,kk) = bin
        cell_vol_sph(mm,bin)  = cell_vol_sph(mm,bin) + volnp(kk)
        dist_from_np(mm,kk)   = r_center_surf
    enddo
enddo
!-----------------------------------------------------------------------------------------------------------!
end subroutine alloc_hist
