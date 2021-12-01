subroutine init_time(ns, ds_ave, ds, koeff)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns

real(8), intent(in)                   :: ds_ave
real(8), intent(out), dimension(ns+1) :: ds, koeff
!------------------------------------------------------------------------------------------------------!
!constant step size scheme
ds = ds_ave

!call simpsonkoef_s(ds_ave, ns, koeff)
call quadinterp_koef(ds, ns, koeff)
!------------------------------------------------------------------------------------------------------!
end subroutine init_time
