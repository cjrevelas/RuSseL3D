subroutine init_arrays()
!-----------------------------------------------------------------------------------------------------------!
use geometry
use kcw
use parser_vars
use arrays
use delta
use error_handing
use iofiles
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
allocate(rdiag1(numnp))
allocate(d2phi_dr2(numnp))
allocate(wa(numnp),wa_mix(numnp),wa_new(numnp),Ufield(numnp))
allocate(ds_matrix_ed(ns_matrix_ed+1),ds_matrix_conv(ns_matrix_conv+1))
allocate(xs_matrix_ed(ns_matrix_ed+1),xs_matrix_conv(ns_matrix_conv+1))
allocate(koeff_matrix_ed(ns_matrix_ed+1),koeff_matrix_conv(ns_matrix_conv+1))
allocate(phia_mx(numnp),phia_gr(numnp),phi_total(numnp))
allocate(qm(numnp,2))
allocate(qm_final(numnp,ns_matrix_ed+1))
allocate(qm_interp_mm(numnp,ns_matrix_conv+1))
allocate(qm_interp_mg(numnp,ns_gr_conv+1))
allocate(qgr_final(numnp,ns_gr_ed+1))
allocate(qgr_interp(numnp,ns_gr_conv+1))

d2phi_dr2         = 0.d0
wa                = 0.d0
wa_mix            = 0.d0
wa_new            = 0.d0
Ufield            = 0.d0
rdiag1            = 0.d0
ds_matrix_ed      = 0.d0
ds_matrix_conv    = 0.d0
xs_matrix_ed      = 0.d0
xs_matrix_conv    = 0.d0
koeff_matrix_ed   = 0.d0
koeff_matrix_conv = 0.d0
phia_mx           = 0.d0
phia_gr           = 0.d0
phi_total         = 0.d0
qm_final          = 0.d0
qm                = 0.d0
qm_interp_mm      = 0.d0
qm_interp_mg      = 0.d0
qgr_final         = 0.d0
qgr_interp        = 0.d0

if (grafted_exist.eq.1) then
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
!-----------------------------------------------------------------------------------------------------------!
end subroutine init_arrays
