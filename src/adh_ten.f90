subroutine adhesion_tension(qf_final, qgr_final, wa, Ufield, phia_fr, phia_gr, part_func, adh_ten)
!-------------------------------------------------------------------------------------------------!
use parser_vars
use geometry
use constants
!-------------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------------!
integer :: k1

real(8), intent(in), dimension(numnp)               :: wa, Ufield, phia_fr, phia_gr
real(8), intent(in), dimension(numnp,ns_free_max+1) :: qf_final
real(8), intent(in), dimension(numnp,ns_gr+1)       :: qgr_final
real(8), intent(in)                                 :: part_func
real(8), intent(out)                                :: adh_ten
real(8)                                             :: Q, sum_final
real(8)                                             :: term1, term2, term3, term4
real(8)                                             :: part_sum1, part_sum2, part_sum3, part_sum4
real(8), dimension(numnp)                           :: dterm1, dterm2
!-------------------------------------------------------------------------------------------------!
term1 = 0.d00
term2 = 0.d00
term3 = 0.d00
term4 = 0.d00

do k1 = 1, numnp
   dterm1(k1) = 0.5d0 * kapa * ((1.d0 - phia_fr(k1) - phia_gr(k1))**2.d0) * chainlen_free !blue: particle-particle interaction
   dterm2(k1) = -(wa(k1) - Ufield(k1)) * (phia_fr(k1) + phia_gr(k1)) * chainlen_free      !pink: rho-w interaction
enddo

call spat_3d(dterm1, term1, Q)
call spat_3d(dterm2, term2, Q)

term1 = term1 * 1.0d-30
term2 = term2 * 1.0d-30
term3 = volume * 1.0d-30 * (1.d00 - part_func) !red: translational entropy of free chains

term4 = 0.d0
do k1 = 1, numnp
    if (qgr_final(k1,1).gt.tol) then
        term4 = term4 + log(qf_final(k1,ns_gr+1))
    endif
enddo
term4 = -chainlen_free / rho_0 * term4 !green: translational entropy of grafted chains

! x1.0D+03 ----> N/m --> mN/m
part_sum1 = term1 * rho_0 * boltz_const_Joule_molK * Temp / chainlen_free * 1.0D+03
part_sum2 = term2 * rho_0 * boltz_const_Joule_molK * Temp / chainlen_free * 1.0D+03
part_sum3 = term3 * rho_0 * boltz_const_Joule_molK * Temp / chainlen_free * 1.0D+03
part_sum4 = term4 * rho_0 * boltz_const_Joule_K * Temp / chainlen_free * 1.0D+03

sum_final = part_sum1 + part_sum2 + part_sum3 + part_sum4

adh_ten = - sum_final/(interf_area*1.D-20) ! adh_ten = (Omega_bulk - Omega) / A

!flush output
open(unit=837, file = 'energy_terms.out.txt', position = 'append')
write(837,'(6(E19.9e2,1X))')  part_sum1, part_sum2, part_sum3, part_sum4, &
   &                          sum_final, adh_ten
close(837)

return
!-------------------------------------------------------------------------------------------------!
end subroutine adhesion_tension
