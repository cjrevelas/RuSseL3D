subroutine get_nchains(numnp, chainlen, rho_mol_bulk, phia, nch_gr)
!-------------------------------------------------------------------------------------------!
use constants, only : n_avog
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

call spat_3d(phia, sum_f, QQ, vol)

nch_gr = sum_f * 1.0d-30 * rho_mol_bulk * n_avog / chainlen

return
!-------------------------------------------------------------------------------------------!
end subroutine get_nchains   
