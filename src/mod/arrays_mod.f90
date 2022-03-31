!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module arrays_mod
!----------------------------------------------------------------------------------------------------------------------------!
real(8), allocatable, dimension(:)   :: ds_mx_ed, ds_gr_ed, ds_mx_conv, ds_gr_conv
real(8), allocatable, dimension(:)   :: coeff_mx_ed, coeff_gr_ed, coeff_mx_conv, coeff_gr_conv
real(8), allocatable, dimension(:)   :: xs_mx_ed, xs_mx_conv, xs_gr_ed, xs_gr_conv
real(8), allocatable, dimension(:)   :: ww, ww_new, ww_mix, Ufield
real(8), allocatable, dimension(:)   :: phi_mx, phi_gr, phi_total
real(8), allocatable, dimension(:)   :: volnp
real(8), allocatable, dimension(:,:) :: phi_gr_indiv
real(8), allocatable, dimension(:)   :: d2phi_dr2
real(8), allocatable, dimension(:,:) :: qmx, qmx_final, qmx_interp_mm, qmx_interp_mg
real(8), allocatable, dimension(:,:) :: qgr, qgr_final, qgr_interp
!----------------------------------------------------------------------------------------------------------------------------!
end module arrays_mod
