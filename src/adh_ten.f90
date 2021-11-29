subroutine adhesion_tension(q_final, wa, Ufield, phia_fr, phia_gr, part_func, adh_ten)
!---------------------------------------------------------------------------------------!
use xdata
use constants
!---------------------------------------------------------------------------------------!
implicit none
!---------------------------------------------------------------------------------------!
integer :: k1

real(8), intent(in), dimension(numnp)      :: wa, Ufield, phia_fr, phia_gr
real(8), intent(in), dimension(numnp,ns+1) :: q_final
real(8), intent(in)                        :: part_func
real(8), intent(out)                       :: adh_ten
real(8)                                    :: Q, L1, L2, sum_final
real(8)                                    :: term1, term2, term3, term4
real(8)                                    :: part_sum1, part_sum2, part_sum3, part_sum4
real(8), dimension(numnp)                  :: dterm1, dterm2
!---------------------------------------------------------------------------------------!
term1 = 0.d00
term2 = 0.d00
term3 = 0.d00
term4 = 0.d00

do k1 = 1, numnp
   dterm1(k1) = 0.5d0 * kapa * ((1.d0 - phia_fr(k1) - phia_gr(k1))**2.d0) !blue: particle-particle interaction
   dterm2(k1) = -(wa(k1) - Ufield(k1)) * (phia_fr(k1) + phia_gr(k1))      !pink: rho-w interaction
enddo

call spat_3d(dterm1, term1, Q)
call spat_3d(dterm2, term2, Q)

term1 = term1 * 1.0d-30
term2 = term2 * 1.0d-30

term3 = volume * 1.0d-30 * (1.d00 - part_func)     !red:   translational entropy of free chains

term4 = -chainlen / rho_0 * log(q_final(1,ns+1))   !green: translational entropy of grafted chains

! x1.0D+03 ----> N/m --> mN/m
part_sum1 = term1 * rho_0 * boltz_const_Joule_molK * Temp / chainlen * 1.0D+03
part_sum2 = term2 * rho_0 * boltz_const_Joule_molK * Temp / chainlen * 1.0D+03
part_sum3 = term3 * rho_0 * boltz_const_Joule_molK * Temp / chainlen * 1.0D+03
part_sum4 = term4 * rho_0 * boltz_const_Joule_molK * Temp / chainlen * 1.0D+03

sum_final   = part_sum1 + part_sum2 + part_sum3! + part_sum4

L1 = box_len(mod(prof_dim+0,3)+1)
L2 = box_len(mod(prof_dim+1,3)+1) 
adh_ten = - sum_final/(2*L1*L2*1.D-20)   ! adh_ten = (Omega_bulk - Omega) / A

return
!---------------------------------------------------------------------------------------!
end subroutine adhesion_tension
