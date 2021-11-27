subroutine adhesion_tension
!--------------------------------------------------------------------!  
use xdata
use constants
!--------------------------------------------------------------------!  
implicit none
!--------------------------------------------------------------------!
real(8)                   :: Q
real(8), dimension(numnp) :: dterm1, dterm2
!--------------------------------------------------------------------!
term1 = 0.d00
term2 = 0.d00
term3 = 0.d00
term4 = 0.d00

do k1 = 1, numnp
   dterm1(k1) = 0.5d0 * kapa * ((1.d0 - phia_new(k1))**2.d0) !blue: particle-particle interaction
   dterm2(k1) = -(wa(k1) - Ufield(k1)) * phia_new(k1)        !pink: rho-w interaction
enddo

call spat_3d(dterm1, term1, Q)
call spat_3d(dterm2, term2, Q)

term1 = term1 * 1.0d-30
term2 = term2 * 1.0d-30
!APS BUGFIX: xc(1,numnp) seems suspecius! Check this!
term3 = volume * 1.0d-30 * (1.d00 - part_func)      !red:   translational entropy of free chains

term4 = -chainlen / rho_0 * log(qf_final(1,ns+1))    !green: translational entropy of grafted chains

! x1.0D+03 ----> N/m --> mN/m
part_sum1 = term1 * rho_0 * boltz_const_Joule_molK * Temp / chainlen * 1.0D+03
part_sum2 = term2 * rho_0 * boltz_const_Joule_molK * Temp / chainlen * 1.0D+03
part_sum3 = term3 * rho_0 * boltz_const_Joule_molK * Temp / chainlen * 1.0D+03
part_sum4 = term4 * rho_0 * boltz_const_Joule_molK * Temp / chainlen * 1.0D+03

sum_final   = part_sum1 + part_sum2 + part_sum3! + part_sum4
adh_ten = - sum_final   ! adh_test = (Omega_bulk - Omega) / A

return
!--------------------------------------------------------------------!  
end subroutine adhesion_tension
