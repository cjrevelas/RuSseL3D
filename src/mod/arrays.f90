module arrays
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
real(8), allocatable, dimension(:)   :: ds_mx_ed, ds_gr_ed, ds_mx_conv, ds_gr_conv
real(8), allocatable, dimension(:)   :: coeff_mx_ed, coeff_gr_ed, coeff_mx_conv, coeff_gr_conv
real(8), allocatable, dimension(:)   :: xs_mx_ed, xs_mx_conv, xs_gr_ed, xs_gr_conv
real(8), allocatable, dimension(:)   :: wa, wa_new, wa_mix, Ufield
real(8), allocatable, dimension(:)   :: phia_mx, phia_gr, phi_total
real(8), allocatable, dimension(:)   :: volnp
real(8), allocatable, dimension(:,:) :: phia_gr_indiv
real(8), allocatable, dimension(:)   :: d2phi_dr2
real(8), allocatable, dimension(:,:) :: qmx, qmx_final, qmx_interp_mm, qmx_interp_mg
real(8), allocatable, dimension(:,:) :: qgr, qgr_final, qgr_interp
!-----------------------------------------------------------------------------------------------------------!
end module arrays
