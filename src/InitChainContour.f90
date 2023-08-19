!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine InitChainContour(contourDiscretization, chainLength, xs_crit, ns, ds_ave, ds, xs, coeff)
!------------------------------------------------------------------------------------------------------!
use constants_mod, only: pi
use flags_mod,     only: contourUniform, contourSymmetric, contourAsymmetric, contourHybrid
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: contourDiscretization, ns
integer             :: kk, ns_part1

real(8), intent(in)                   :: chainLength, ds_ave, xs_crit
real(8), intent(out), dimension(ns+1) :: ds, xs, coeff
!------------------------------------------------------------------------------------------------------!
if (contourDiscretization.eq.contourUniform) then
  ds = ds_ave
  do kk = 2, ns+1
    xs(kk) = xs(kk-1) + ds_ave
  enddo
elseif (contourDiscretization.eq.contourSymmetric) then
  ds(1)=0.0d0
  do kk = 2, ns+1
    xs(kk) = chainLength * 0.5d0 * (1.0d0 - DCOS(pi * (DBLE(kk)-1.0d0) /  DBLE(ns)))
    ds(kk) = xs(kk) - xs(kk-1)
  enddo
elseif (contourDiscretization.eq.contourAsymmetric) then
  ds(1)=0.0d0
  do kk = 2, ns+1
    xs(kk) = chainLength *         (1.0d0 - DCOS(pi * (DBLE(kk)-1.0d0) / (DBLE(ns) * 2.0d0)))
    ds(kk) = xs(kk) - xs(kk-1)
  enddo
elseif (contourDiscretization.eq.contourHybrid) then
  ds(1)=0.0d0
  ns_part1 = 2 * NINT(0.5d0 * xs_crit / ds_ave)
  if (xs_crit.le.chainLength) then
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

call ComputeContourCoeffs(ds, ns, coeff)
return
!------------------------------------------------------------------------------------------------------!
end subroutine InitChainContour
