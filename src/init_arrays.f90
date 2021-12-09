subroutine init_arrays()
!-----------------------------------------------------------------------------------------------------------!
use geometry,    only: numnp
use kcw,         only: rdiag1
use parser_vars, only: ns_mx_ed, ns_mx_conv, mx_exist, gr_exist, ns_gr_ed, ns_gr_conv
use arrays
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
allocate(qm(2,numnp))
allocate(qm_final(ns_matrix_ed+1,numnp))
allocate(qm_interp_mm(ns_matrix_conv+1,numnp))
allocate(qm_interp_mg(ns_gr_conv+1,numnp))
allocate(qgr_final(ns_gr_ed+1,numnp))
allocate(qgr_interp(ns_gr_conv+1,numnp))

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
    allocate(qgr(2,numnp))

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
