subroutine init_chain_contour(sym, chainlen, ns, ds_ave, ds, xs, koeff)
!------------------------------------------------------------------------------------------------------!
use constants
use parser_vars
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns

logical, intent(in) :: sym

real(8), intent(in)                   :: chainlen, ds_ave
real(8), intent(out), dimension(ns+1) :: ds, xs, koeff

integer :: k1
!------------------------------------------------------------------------------------------------------!
!constant step size scheme
if (contour_integration_scheme.eq.1) then
    ds = ds_ave
    do k1 = 2, ns+1
        xs(k1) = xs(k1-1) + ds_ave
    enddo
elseif (contour_integration_scheme.eq.2) then
    ds(1)=0.d0
    do k1 = 2, ns+1
        if (sym) then
            xs(k1) = chainlen * 0.5d0 * (1.d0 - DCOS(pi * (dble(k1)-1.d0) /  dble(ns)))         !symmetric scheme
        else
            xs(k1) = chainlen *         (1.d0 - DCOS(pi * (dble(k1)-1.d0) / (dble(ns) * 2.d0))) !asymmetric scheme
        endif

    ds(k1) = xs(k1) - xs(k1-1)
    enddo
endif

call quadinterp_koef(ds, ns, koeff)
!------------------------------------------------------------------------------------------------------!
end subroutine init_chain_contour
