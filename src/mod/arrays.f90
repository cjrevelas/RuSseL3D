module arrays
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
real(8), allocatable, dimension(:)   :: ds_matrix_ed, ds_gr_ed, ds_matrix_conv, ds_gr_conv
real(8), allocatable, dimension(:)   :: koeff_matrix_ed, koeff_gr_ed, koeff_matrix_conv, koeff_gr_conv
real(8), allocatable, dimension(:)   :: xs_matrix_ed, xs_matrix_conv, xs_gr_ed, xs_gr_conv
real(8), allocatable, dimension(:)   :: wa, wa_new, wa_mix, Ufield
real(8), allocatable, dimension(:)   :: phia_mx, phia_gr, phi_total
real(8), allocatable, dimension(:,:) :: qm, qm_final, qm_interp_mm, qm_interp_mg
real(8), allocatable, dimension(:,:) :: qgr, qgr_final, qgr_interp
!-----------------------------------------------------------------------------------------------------------!
end module arrays
