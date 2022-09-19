!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ComputeNumberOfChains(chainLength, phi, numChains)
!-------------------------------------------------------------------------------------------!
use geometry_mod,    only: numNodes
use parser_vars_mod, only: molarBulkDensity
use constants_mod,   only: n_avog, A3_to_m3
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
real(8), intent(in)                      :: chainLength
real(8), intent(in), dimension(numNodes) :: phi
real(8), intent(out)                     :: numChains
real(8)                                  :: sum_f, QQ, vol
!-------------------------------------------------------------------------------------------!
sum_f = 0.0d0

call FemIntegration(phi, sum_f, QQ, vol)

numChains = sum_f * A3_to_m3 * molarBulkDensity * n_avog / chainLength

return
!-------------------------------------------------------------------------------------------!
end subroutine ComputeNumberOfChains
