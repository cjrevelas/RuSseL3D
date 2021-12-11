!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

real(8) function compute_stretching_energy(gnode_id, qmx, qgr)
!-------------------------------------------------------------------------------------------------!
use constants_mod,   only: A3_to_m3
use parser_vars_mod, only: beta, ns_gr_conv, chainlen_gr, Rg2_per_mon_gr, rho_seg_bulk
use geometry_mod,    only: numnp, xc
!-------------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------------!
integer, intent(in) :: gnode_id
integer             :: ii

real(8), intent(in), dimension(ns_gr_conv+1,numnp) :: qmx, qgr
real(8), dimension(numnp)                          :: dr2, phi_end, rho_end, A_stretch
real(8)                                            :: Q, vol
!-------------------------------------------------------------------------------------------------!
do ii = 1, numnp
    phi_end(ii) = 1 / chainlen_gr * qgr(ns_gr_conv+1,ii) * qmx(1,ii)
    rho_end(ii) = phi_end(ii) * rho_seg_bulk * A3_to_m3
    dr2(ii)      = (xc(1,ii)-xc(1,gnode_id))**2 + (xc(2,ii)-xc(2,gnode_id))**2 + (xc(3,ii)-xc(3,gnode_id))**2

    A_stretch(ii) = 3.0/(2.d0 * beta * (Rg2_per_mon_gr*chainlen_gr*6.d0)) * dr2(ii)
enddo
call fem_integration(rho_end*A_stretch, compute_stretching_energy, Q, vol)
!-------------------------------------------------------------------------------------------------!
end function compute_stretching_energy
