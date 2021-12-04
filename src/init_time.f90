subroutine init_time(ns, ds_ave, ds, xs, koeff)
!------------------------------------------------------------------------------------------------------!
use constants
use parser_vars
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns

real(8), intent(in)                   :: ds_ave
real(8), intent(out), dimension(ns+1) :: ds, xs, koeff

integer :: k1
!------------------------------------------------------------------------------------------------------!
!constant step size scheme
if (time_integration_scheme.eq.1) then
    ds = ds_ave
    do k1 = 2, ns+1
        xs(k1) = xs(k1-1) + ds_ave
    enddo
elseif (time_integration_scheme.eq.2) then
    ds(1)=0.d0
    do k1 = 2, ns+1
!       xs(k1) = 0.5d0 * (1.d0 - DCOS(pi * (dble(k1)-1.d0) /  dble(ns)))         ! Symmetric scheme
        xs(k1) =         (1.d0 - DCOS(pi * (dble(k1)-1.d0) / (dble(ns) * 2.d0))) ! Asymmetric scheme
        ds(k1) = xs(k1) - xs(k1-1)
    enddo
endif

call quadinterp_koef(ds, ns, koeff)
!------------------------------------------------------------------------------------------------------!
end subroutine init_time
