subroutine grafted_chains(numnp, chainlen, rho_0, phia, nch_gr)
!-------------------------------------------------------------------------------------------!
use constants
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in)                   :: numnp

real(8), intent(in)                   :: chainlen, rho_0
real(8), intent(in), dimension(numnp) :: phia
real(8), intent(out)                  :: nch_gr
real(8)                               :: sum_f, Q
!-------------------------------------------------------------------------------------------!
sum_f = 0.d0

call spat_3d(phia, sum_f, Q)

nch_gr = sum_f * 1.0d-30 * rho_0 * avogadro_constant / chainlen

return
!-------------------------------------------------------------------------------------------!
end subroutine grafted_chains
