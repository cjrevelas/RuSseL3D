!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_arrays()
!-----------------------------------------------------------------------------------------------------------!
use geometry_mod,    only: numNodes
use kcw_mod,         only: rdiag1
use parser_vars_mod, only: numEdwPointsMatrix, numConvolPointsMatrix, matrixExist, graftedExist, &
                           numEdwPointsGrafted, numConvolPointsGrafted
use arrays_mod
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
allocate(rdiag1(numNodes))
allocate(d2phi_dr2(numNodes))
allocate(ww(numNodes),ww_mix(numNodes),ww_new(numNodes),Ufield(numNodes))
allocate(phi_mx(numNodes),phi_gr(numNodes),phi_total(numNodes))
allocate(volnp(numNodes))
allocate(ds_mx_ed(numEdwPointsMatrix+1))
allocate(xs_mx_ed(numEdwPointsMatrix+1))
allocate(coeff_mx_ed(numEdwPointsMatrix+1))
allocate(qmx(2,numNodes))
allocate(qmx_final(numEdwPointsMatrix+1,numNodes))
allocate(qmx_interp_mg(numConvolPointsGrafted+1,numNodes))
allocate(qgr_interp(numConvolPointsGrafted+1,numNodes))

if (matrixExist.eq.1) then
  allocate(qmx_interp_mm(numConvolPointsMatrix+1,numNodes))
  allocate(ds_mx_conv(numConvolPointsMatrix+1))
  allocate(xs_mx_conv(numConvolPointsMatrix+1))
  allocate(coeff_mx_conv(numConvolPointsMatrix+1))

  qmx_interp_mm = 0.0d0
  ds_mx_conv    = 0.0d0
  xs_mx_conv    = 0.0d0
  coeff_mx_conv = 0.0d0
endif

volnp         = 0.0d0
d2phi_dr2     = 0.0d0
ww            = 0.0d0
ww_mix        = 0.0d0
ww_new        = 0.0d0
Ufield        = 0.0d0
rdiag1        = 0.0d0
phi_mx        = 0.0d0
phi_gr        = 0.0d0
phi_total     = 0.0d0
ds_mx_ed      = 0.0d0
xs_mx_ed      = 0.0d0
coeff_mx_ed   = 0.0d0
qmx           = 0.0d0
qmx_final     = 0.0d0
qmx_interp_mg = 0.0d0

if (graftedExist.eq.1) then
  allocate(ds_gr_ed(numEdwPointsGrafted+1),ds_gr_conv(numConvolPointsGrafted+1))
  allocate(xs_gr_ed(numEdwPointsGrafted+1),xs_gr_conv(numConvolPointsGrafted+1))
  allocate(coeff_gr_ed(numEdwPointsGrafted+1),coeff_gr_conv(numConvolPointsGrafted+1))
  allocate(qgr(2,numNodes))
  allocate(qgr_final(numEdwPointsGrafted+1,numNodes))

  ds_gr_ed      = 0.0d0
  ds_gr_conv    = 0.0d0
  xs_gr_ed      = 0.0d0
  xs_gr_conv    = 0.0d0
  coeff_gr_ed   = 0.0d0
  coeff_gr_conv = 0.0d0
  qgr           = 0.0d0
  qgr_final     = 0.0d0
  qgr_interp    = 0.0d0
endif
!-----------------------------------------------------------------------------------------------------------!
end subroutine init_arrays
