subroutine init_time
!------------------------------------------------------------------------------------------------------!
use xdata
use constants
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: k1
!------------------------------------------------------------------------------------------------------!
!discretize time domain
allocate(ds(ns+1))
allocate(koeff(ns+1))
allocate(xs(ns+1))
ds    = 0.d0
koeff = 0.d0
xs    = 0.d0

!constant step size scheme
if (time_integration_scheme.eq.1) then
    ds = ds_ave
    call simpsonkoef_s(ds_ave, ns, koeff)
endif

!non constant step size scheme with quadradic interpolation
if (time_integration_scheme.eq.2) then

! A)
!    ds = ds_ave
!    do k1 = 2, ns+1
!        xs(k1) = xs(k1-1) + ds_ave
!        ds(k1) = xs(k1) - xs(k1-1)
!    enddo
! B)
    ds(1)=0.d0
    do k1 = 2, ns+1
        xs(k1) = 0.5d0 * (1.d0 - DCOS(pi * (dble(k1)-1.d0) / dble(ns)))
        ds(k1) = xs(k1) - xs(k1-1)
    enddo

    call quadinterp_koef(ds, ns, xs, koeff)
endif
!------------------------------------------------------------------------------------------------------!
end subroutine init_time
