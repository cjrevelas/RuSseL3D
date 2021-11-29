subroutine grafted_chains(phia, nch_per_area)
!-------------------------------------------------------------------------------------------!
use xdata
use constants
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
real(8), intent(in), dimension(numnp)      :: phia
real(8), intent(out)                       :: nch_per_area
real(8)                                    :: sum_f, Q
!-------------------------------------------------------------------------------------------!
sum_f = 0.d0

call spat_3d(phia, sum_f, Q)

nch_per_area = sum_f*1.0d-30*rho_0*avogadro_constant/chainlen

return
!-------------------------------------------------------------------------------------------!
end subroutine grafted_chains
