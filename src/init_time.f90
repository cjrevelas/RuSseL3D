subroutine init_time(ds, koeff)
!------------------------------------------------------------------------------------------------------!
use xdata
use constants
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(out), dimension(ns+1) :: ds, koeff
!------------------------------------------------------------------------------------------------------!
!constant step size scheme
ds = ds_ave
call simpsonkoef_s(ds_ave, ns, koeff)
!------------------------------------------------------------------------------------------------------!
end subroutine init_time
