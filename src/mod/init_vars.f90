module init_vars
!-----------------------------------------------------------------------------------------------------------!
use geometry
use kcw
use parser_vars
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
real(8), allocatable, dimension(:)   :: ds_matrix_ed, ds_gr_ed, ds_matrix_conv, ds_gr_conv
real(8), allocatable, dimension(:)   :: koeff_matrix_ed, koeff_gr_ed, koeff_matrix_conv, koeff_gr_conv
real(8), allocatable, dimension(:)   :: xs_matrix_ed, xs_matrix_conv, xs_gr_ed, xs_gr_conv
real(8), allocatable, dimension(:)   :: wa, wa_new, wa_mix, Ufield
real(8), allocatable, dimension(:)   :: phia_mx, phia_gr
real(8), allocatable, dimension(:,:) :: qm, qm_final, qm_interp_mm, qm_interp_mg
real(8), allocatable, dimension(:,:) :: qgr, qgr_final, qgr_interp
!-----------------------------------------------------------------------------------------------------------!
contains

subroutine init_vars_and_arrays
implicit none

allocate(rdiag1(numnp))
allocate(wa(numnp),wa_mix(numnp),wa_new(numnp),Ufield(numnp))
allocate(ds_matrix_ed(ns_matrix_ed+1),ds_matrix_conv(ns_matrix_conv+1))
allocate(xs_matrix_ed(ns_matrix_ed+1),xs_matrix_conv(ns_matrix_conv+1))
allocate(koeff_matrix_ed(ns_matrix_ed+1),koeff_matrix_conv(ns_matrix_conv+1))
allocate(phia_mx(numnp))
allocate(qm(numnp,2))
allocate(qm_final(numnp,ns_matrix_ed+1))
allocate(qm_interp_mm(numnp,ns_matrix_conv+1))
allocate(qm_interp_mg(numnp,ns_gr_conv+1))
allocate(qgr_final(numnp,ns_gr_ed+1))
allocate(qgr_interp(numnp,ns_gr_conv+1))
allocate(phia_gr(numnp))

wa     = 0.d0
wa_mix = 0.d0
wa_new = 0.d0
Ufield = 0.d0
rdiag1 = 0.d0
ds_matrix_ed      = 0.d0
ds_matrix_conv    = 0.d0
xs_matrix_ed      = 0.d0
xs_matrix_conv    = 0.d0
koeff_matrix_ed   = 0.d0
koeff_matrix_conv = 0.d0
phia_mx           = 0.d0
qm_final          = 0.d0
qm                = 0.d0
qm_interp_mm      = 0.d0
qm_interp_mg      = 0.d0
qgr_final  = 0.d0
qgr_interp = 0.d0
phia_gr = 0.d0

if (use_grafted.eq.1) then
    allocate(ds_gr_ed(ns_gr_ed+1),ds_gr_conv(ns_gr_conv+1))
    allocate(xs_gr_ed(ns_gr_ed+1),xs_gr_conv(ns_gr_conv+1))
    allocate(koeff_gr_ed(ns_gr_ed+1),koeff_gr_conv(ns_gr_conv+1))
    allocate(qgr(numnp,2))

    ds_gr_ed      = 0.d0
    ds_gr_conv    = 0.d0
    xs_gr_ed      = 0.d0
    xs_gr_conv    = 0.d0
    koeff_gr_ed   = 0.d0
    koeff_gr_conv = 0.d0
    qgr           = 0.d0
endif

end subroutine init_vars_and_arrays
!-----------------------------------------------------------------------------------------------------------!
end module init_vars
