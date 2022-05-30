!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_chain_contour(contour_discr, chainlen, xs_crit, ns, ds_ave, ds, xs, coeff)
!------------------------------------------------------------------------------------------------------!
use constants_mod, only: pi
use flags_mod,     only: contour_uniform, contour_symm, contour_asymm, contour_hybrid
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: contour_discr, ns
integer             :: kk, ns_part1

real(8), intent(in)                   :: chainlen, ds_ave, xs_crit
real(8), intent(out), dimension(ns+1) :: ds, xs, coeff
!------------------------------------------------------------------------------------------------------!
if (contour_discr.eq.contour_uniform) then
  ds = ds_ave
  do kk = 2, ns+1
    xs(kk) = xs(kk-1) + ds_ave
  enddo
elseif (contour_discr.eq.contour_symm) then
  ds(1)=0.0d0
  do kk = 2, ns+1
    xs(kk) = chainlen * 0.5d0 * (1.0d0 - DCOS(pi * (DBLE(kk)-1.0d0) /  DBLE(ns)))
    ds(kk) = xs(kk) - xs(kk-1)
  enddo
elseif (contour_discr.eq.contour_asymm) then
  ds(1)=0.0d0
  do kk = 2, ns+1
    xs(kk) = chainlen *         (1.0d0 - DCOS(pi * (DBLE(kk)-1.0d0) / (DBLE(ns) * 2.0d0)))
    ds(kk) = xs(kk) - xs(kk-1)
  enddo
elseif (contour_discr.eq.contour_hybrid) then
  ds(1)=0.0d0
  ns_part1 = 2 * NINT(0.5d0 * xs_crit / ds_ave)
  if (xs_crit.le.chainlen) then
    do kk = 2, ns_part1+1
      xs(kk) = xs_crit * (1.0d0 - DCOS(pi * (DBLE(kk)-1.0d0) / (DBLE(ns_part1) * 2.0d0)))
      ds(kk) = xs(kk)  - xs(kk-1)
    enddo

    do kk = ns_part1+2, ns+1
      xs(kk) = xs(kk-1) + ds(ns_part1+1)
      ds(kk) = xs(kk) - xs(kk-1)
    enddo
  else
    do kk = 2, ns + 1
      xs(kk) = xs_crit * (1.0d0 - DCOS(pi * (DBLE(kk)-1.0d0) / (DBLE(ns_part1) * 2.0d0)))
      ds(kk) = xs(kk)  - xs(kk-1)
    enddo
  endif
endif

call compute_contour_coeffs(ds, ns, coeff)
return
!------------------------------------------------------------------------------------------------------!
end subroutine init_chain_contour
