!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine InitChainContour(contourDiscretization, chainLength, xsCrit, ns, dsAve, ds, xs, coeff)
!------------------------------------------------------------------------------------------------------!
use constants_mod, only: pi
use flags_mod,     only: contourUniform, contourSymmetric, contourAsymmetric, contourHybrid
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: contourDiscretization, ns
integer             :: kk, nsPartOne

real(8), intent(in)                   :: chainLength, dsAve, xsCrit
real(8), intent(out), dimension(ns+1) :: ds, xs, coeff
!------------------------------------------------------------------------------------------------------!
if (contourDiscretization.eq.contourUniform) then
  ds = dsAve
  do kk = 2, ns+1
    xs(kk) = xs(kk-1) + dsAve
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
    xs(kk) = chainLength * (1.0d0 - DCOS(pi * (DBLE(kk)-1.0d0) / (DBLE(ns) * 2.0d0)))
    ds(kk) = xs(kk) - xs(kk-1)
  enddo
elseif (contourDiscretization.eq.contourHybrid) then
  ds(1)=0.0d0
  nsPartOne = 2 * NINT(0.5d0 * xsCrit / dsAve)

  if (xsCrit.le.chainLength) then
    do kk = 2, nsPartOne+1
      xs(kk) = xsCrit * (1.0d0 - DCOS(pi * (DBLE(kk)-1.0d0) / (DBLE(nsPartOne) * 2.0d0)))
      ds(kk) = xs(kk)  - xs(kk-1)
    enddo

    do kk = nsPartOne+2, ns+1
      xs(kk) = xs(kk-1) + ds(nsPartOne+1)
      ds(kk) = xs(kk) - xs(kk-1)
    enddo
  else
    do kk = 2, ns + 1
      xs(kk) = xsCrit * (1.0d0 - DCOS(pi * (DBLE(kk)-1.0d0) / (DBLE(nsPartOne) * 2.0d0)))
      ds(kk) = xs(kk)  - xs(kk-1)
    enddo
  endif
endif

CALL ComputeContourCoeffs(ds, ns, coeff)

return
!------------------------------------------------------------------------------------------------------!
end subroutine InitChainContour
