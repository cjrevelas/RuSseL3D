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
allocate(phia_mx(numnp),phia_gr(numnp),phi_total(numnp))
allocate(volnp(numnp))
allocate(ds_mx_ed(ns_mx_ed+1))
allocate(xs_mx_ed(ns_mx_ed+1))
allocate(coeff_mx_ed(ns_mx_ed+1))
allocate(qmx(2,numnp))
allocate(qmx_final(ns_mx_ed+1,numnp))
allocate(qmx_interp_mg(ns_gr_conv+1,numnp))
allocate(qgr_interp(ns_gr_conv+1,numnp))

if (mx_exist.eq.1) then
    allocate(qmx_interp_mm(ns_mx_conv+1,numnp))
    allocate(ds_mx_conv(ns_mx_conv+1))
    allocate(xs_mx_conv(ns_mx_conv+1))
    allocate(coeff_mx_conv(ns_mx_conv+1))

    qmx_interp_mm = 0.d0
    ds_mx_conv    = 0.d0
    xs_mx_conv    = 0.d0
    coeff_mx_conv = 0.d0
endif

volnp         = 0.d0
d2phi_dr2     = 0.d0
wa            = 0.d0
wa_mix        = 0.d0
wa_new        = 0.d0
Ufield        = 0.d0
rdiag1        = 0.d0
phia_mx       = 0.d0
phia_gr       = 0.d0
phi_total     = 0.d0
ds_mx_ed      = 0.d0
xs_mx_ed      = 0.d0
coeff_mx_ed   = 0.d0
qmx           = 0.d0
qmx_final     = 0.d0
qmx_interp_mg = 0.d0

if (gr_exist.eq.1) then
    allocate(ds_gr_ed(ns_gr_ed+1),ds_gr_conv(ns_gr_conv+1))
    allocate(xs_gr_ed(ns_gr_ed+1),xs_gr_conv(ns_gr_conv+1))
    allocate(coeff_gr_ed(ns_gr_ed+1),coeff_gr_conv(ns_gr_conv+1))
    allocate(qgr(2,numnp))
    allocate(qgr_final(ns_gr_ed+1,numnp))

    ds_gr_ed      = 0.d0
    ds_gr_conv    = 0.d0
    xs_gr_ed      = 0.d0
    xs_gr_conv    = 0.d0
    coeff_gr_ed   = 0.d0
    coeff_gr_conv = 0.d0
    qgr           = 0.d0
    qgr_final     = 0.d0
    qgr_interp    = 0.d0
endif
!-----------------------------------------------------------------------------------------------------------!
end subroutine init_arrays
