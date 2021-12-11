!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_number_of_chains(numnp, chainlen, rho_mol_bulk, phia, num_chains)
!-------------------------------------------------------------------------------------------!
use constants_mod, only : n_avog, A3_to_m3
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in)                   :: numnp

real(8), intent(in)                   :: chainlen, rho_mol_bulk
real(8), intent(in), dimension(numnp) :: phia
real(8), intent(out)                  :: num_chains
real(8)                               :: sum_f, QQ, vol
!-------------------------------------------------------------------------------------------!
sum_f = 0.d0

call fem_integration(phia, sum_f, QQ, vol)

num_chains = sum_f * A3_to_m3 * rho_mol_bulk * n_avog / chainlen

return
!-------------------------------------------------------------------------------------------!
end subroutine compute_number_of_chains
