subroutine init_scf_params()
!------------------------------------------------------------------------------------------------------!
use parser_vars,  only: chainlen_mx, iow, k_gr, k_gr_tilde, massden, mon_mass, pres, &
&                       rho_mass_bulk, rho_mol_bulk, rho_seg_bulk, temp
use constants,    only: boltz_const_Joule_K, boltz_const_Joule_molK, gr_cm3_to_kg_m3, N_avog
use eos,          only: eos_type, kapa, V_star, P_star, T_star, rho_star, rho_tilde_bulk, P_tilde, T_tilde, &
&                       rsl_N, hlf_kappa_T, eos_rho_tilde_0
use flags,        only: eos_helfand, eos_sl
use write_helper, only: adjl
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
if (eos_type.eq.eos_helfand) then
    rho_mol_bulk = massden/mon_mass*1.d06
    write(iow,'(3X,A40,E16.9,A10)')adjl("Segment density in bulk (rho):",40),rho_mol_bulk," [mol/m^3]"
    write(6  ,'(3X,A40,E16.9,A10)')adjl("Segment density in bulk (rho):",40),rho_mol_bulk," [mol/m^3]"

    kapa = 1.d0 /(hlf_kappa_T * boltz_const_Joule_molK * Temp * rho_mol_bulk)
    write(iow,'(3X,A40,E16.9)')adjl("kapa = 1/[k_T k_B T rho_0]:",40),kapa
    write(6  ,'(3X,A40,E16.9)')adjl("kapa = 1/[k_T k_B T rho_0]:",40),kapa
elseif (eos_type.eq.eos_sl) then
    write(iow,'(3X,A45)') adjl("Computation of polymer  mass density from SL EoS",45)
    write(*  ,'(3X,A45)') adjl("Computation of polymer  mass density from SL EoS",45)
    V_star         = boltz_const_Joule_K * T_star / P_star
    T_tilde        = Temp / T_star
    P_tilde        = Pres / P_star
    rsl_N          = (mon_mass * P_star) / (rho_star * 1.d03 * boltz_const_Joule_molK * T_star)
    rho_tilde_bulk = eos_rho_tilde_0(T_tilde, P_tilde, rsl_N*chainlen_mx)
    rho_mass_bulk  = rho_tilde_bulk * rho_star
    rho_mol_bulk   = rho_mass_bulk / mon_mass * gr_cm3_to_kg_m3
    rho_seg_bulk   = rho_mol_bulk * N_avog

    write(iow,'(3X,A45,F16.4,'' g/cm3'')')adjl('mass density was recomputed as:',45), rho_mass_bulk/gr_cm3_to_kg_m3
    write(*  ,'(3X,A45,F16.4,'' g/cm3'')')adjl('mass density was recomputed as:',45), rho_mass_bulk/gr_cm3_to_kg_m3

    k_gr = 2.d0 * P_star * rsl_N**2 * V_star**(8.d0/3.d0) * k_gr_tilde
endif
!------------------------------------------------------------------------------------------------------!
end subroutine init_scf_params
