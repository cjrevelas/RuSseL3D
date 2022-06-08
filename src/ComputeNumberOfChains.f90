!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ComputeNumberOfChains(numNodes, chainlen, molarBulkDensity, phi, num_chains)
!-------------------------------------------------------------------------------------------!
use constants_mod, only : n_avog, A3_to_m3
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in) :: numNodes

real(8), intent(in)                      :: chainlen, molarBulkDensity
real(8), intent(in), dimension(numNodes) :: phi
real(8), intent(out)                     :: num_chains
real(8)                                  :: sum_f, QQ, vol
!-------------------------------------------------------------------------------------------!
sum_f = 0.0d0

call FemIntegration(phi, sum_f, QQ, vol)

num_chains = sum_f * A3_to_m3 * molarBulkDensity * n_avog / chainlen

return
!-------------------------------------------------------------------------------------------!
end subroutine ComputeNumberOfChains
