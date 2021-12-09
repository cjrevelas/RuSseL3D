subroutine compute_number_of_chains(numnp, chainlen, rho_mol_bulk, phia, nch_gr)
!-------------------------------------------------------------------------------------------!
use constants, only : n_avog, A3_to_m3
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in)                   :: numnp

real(8), intent(in)                   :: chainlen, rho_mol_bulk
real(8), intent(in), dimension(numnp) :: phia
real(8), intent(out)                  :: nch_gr
real(8)                               :: sum_f, QQ, vol
!-------------------------------------------------------------------------------------------!
sum_f = 0.d0

call fem_integration(phia, sum_f, QQ, vol)

nch_gr = sum_f * A3_to_m3 * rho_mol_bulk * n_avog / chainlen

return
!-------------------------------------------------------------------------------------------!
end subroutine compute_number_of_chains
