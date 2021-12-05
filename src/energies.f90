subroutine energies(qm_interp_mg, qgr_interp, wa, Ufield, phia_mx, phia_gr, part_func, adh_ten)
!-------------------------------------------------------------------------------------------------!
use parser_vars
use geometry
use constants
!-------------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------------!
integer :: k1

real(8), intent(in), dimension(numnp)               :: wa, Ufield, phia_mx, phia_gr
real(8), intent(in), dimension(numnp,ns_gr_conv+1)  :: qm_interp_mg
real(8), intent(in), dimension(numnp,ns_gr_conv+1)  :: qgr_interp
real(8), intent(in)                                 :: part_func
real(8), intent(out)                                :: adh_ten
real(8)                                             :: Q, sum_final, vol
real(8)                                             :: term1, term2, term3, term4
real(8)                                             :: part_sum1, part_sum2, part_sum3, part_sum4
real(8), dimension(numnp)                           :: dterm1, dterm2
!-------------------------------------------------------------------------------------------------!
term1 = 0.d00
term2 = 0.d00
term3 = 0.d00
term4 = 0.d00

do k1 = 1, numnp
   dterm1(k1) = 0.5d0 * kapa * ((1.d0 - phia_mx(k1) - phia_gr(k1))**2.d0) * chainlen_matrix !blue: particle-particle interaction
   dterm2(k1) = -(wa(k1) - Ufield(k1)) * (phia_mx(k1) + phia_gr(k1)) * chainlen_matrix      !pink: rho-w interaction
enddo

call spat_3d(dterm1, term1, Q, vol)
call spat_3d(dterm2, term2, Q, vol)

term1 = term1 * 1.0d-30
term2 = term2 * 1.0d-30
term3 = vol   * 1.0d-30 * (1.d00 - part_func) !red: translational entropy of matrix chains

term4 = 0.d0
do k1 = 1, numnp
    if (qgr_interp(k1,1).gt.tol) then
        term4 = term4 + log(qm_interp_mg(k1,ns_gr_conv+1))
    endif
enddo
term4 = -chainlen_matrix / rho_0 * term4 !green: translational entropy of grafted chains

! x1.0D+03 ----> N/m --> mN/m
part_sum1 = term1 * rho_0 * boltz_const_Joule_molK * Temp / chainlen_matrix * 1.0D+03
part_sum2 = term2 * rho_0 * boltz_const_Joule_molK * Temp / chainlen_matrix * 1.0D+03
part_sum3 = term3 * rho_0 * boltz_const_Joule_molK * Temp / chainlen_matrix * 1.0D+03
part_sum4 = term4 * rho_0 * boltz_const_Joule_K    * Temp / chainlen_matrix * 1.0D+03

sum_final = part_sum1 + part_sum2 + part_sum3 + part_sum4

adh_ten = - sum_final/(interf_area*1.D-20) ! adh_ten = (Omega_bulk - Omega) / A

!flush output
open(unit=837, file = 'energy_terms.out.txt')
write(837,'(6(A19,1X))')      "part_sum1", "part_sum2", "part_sum3", "part_sum4", &
   &                          "sum_final", "adh_ten"
write(837,'(6(E19.9e2,1X))')  part_sum1, part_sum2, part_sum3, part_sum4, &
   &                          sum_final, adh_ten
close(837)

return
!-------------------------------------------------------------------------------------------------!
end subroutine energies
