subroutine energies(qm_interp_mg, qgr_interp, wa, Ufield, phia_mx, phia_gr, part_func, adh_ten)
!-------------------------------------------------------------------------------------------------!
use parser_vars
use geometry
use constants
use iofiles
!-------------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------------!
integer :: k1

real(8), intent(in), dimension(numnp,ns_gr_conv+1)  :: qm_interp_mg, qgr_interp
real(8), intent(in), dimension(numnp)               :: wa, Ufield, phia_mx, phia_gr
real(8), intent(in)                                 :: part_func
real(8), intent(out)                                :: adh_ten
real(8)                                             :: Q=0.d0, vol=0.d0
real(8)                                             :: term1=0.d0, term2=0.d0, term3=0.d0, term4=0.d0
real(8), dimension(numnp)                           :: dterm1, dterm2
!-------------------------------------------------------------------------------------------------!
dterm1 = 0.d0
dterm2 = 0.d0

do k1 = 1, numnp
   dterm1(k1) = (1.d0 - phia_mx(k1) - phia_gr(k1))**2.d0             !blue: particle-particle interaction
   dterm2(k1) = -(wa(k1) - Ufield(k1)) * (phia_mx(k1) + phia_gr(k1)) !pink: rho-w interaction
enddo

call spat_3d(dterm1, term1, Q, vol)
call spat_3d(dterm2, term2, Q, vol)

term1 = term1 * 1.0d-30 * 0.5d0 / kappa_T
term2 = term2 * 1.0d-30 * rho_0 * boltz_const_Joule_molK * Temp
term3 = rho_0 * 1.0d-30 * vol   * boltz_const_Joule_molK * Temp * (1.d00 - part_func) / chainlen_matrix !red: entropy of matrix chains

do k1 = 1, numnp
    if (qgr_interp(k1,1).gt.tol) term4 = term4 + log(qm_interp_mg(k1,ns_gr_conv+1))
enddo
term4 = - boltz_const_Joule_K * Temp * term4 !green: entropy of grafted chains

adh_ten = (term1 + term2 + term3 + term4) * 1.d03 / (interf_area*1.d-20) 

open(unit=837, file = energy_terms)
write(837,'(6(A19,1X))')      "term1", "term2", "term3", "term4", "adh_ten"
write(837,'(6(E19.9e2,1X))')  term1*1.d03, term2*1.d03, term3*1.d03, term4*1.d03, adh_ten
close(837)

return
!-------------------------------------------------------------------------------------------------!
end subroutine energies
