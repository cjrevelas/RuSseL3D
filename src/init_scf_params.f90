!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_scf_params()
!------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: lengthMatrix, iow, sgtParam, sgtParamTilde, massDensity, massOfMonomer, &
                            massBulkDensity, molarBulkDensity, segmentBulkDensity, temperature, pressure
use constants_mod,    only: boltz_const_Joule_K, boltz_const_Joule_molK, gr_cm3_to_kg_m3, N_avog, &
                            kg_m3_to_gr_m3, m3_to_cm3
use eos_mod,          only: eos_type, kapa, V_star, P_star, T_star, rho_star, rho_tilde_bulk, &
                            P_tilde, T_tilde, rsl_N, hlf_kappa_T, eos_rho_tilde_0
use flags_mod,        only: eos_helfand, eos_sl
use write_helper_mod, only: adjl
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
if (eos_type.eq.eos_helfand) then
  molarBulkDensity = massDensity / massOfMonomer * m3_to_cm3
  segmentBulkDensity = molarBulkDensity * N_avog
  write(iow,'(3X,A40,E16.9,A10)')adjl("Segment density in bulk (rho):",40), molarBulkDensity, " [mol/m^3]"
  write(6  ,'(3X,A40,E16.9,A10)')adjl("Segment density in bulk (rho):",40), molarBulkDensity, " [mol/m^3]"

  kapa = 1.0d0 /(hlf_kappa_T * boltz_const_Joule_molK * temperature * molarBulkDensity)
  write(iow,'(3X,A40,E16.9)')adjl("kapa = 1/[k_T k_B T rho_0]:",40), kapa
  write(6  ,'(3X,A40,E16.9)')adjl("kapa = 1/[k_T k_B T rho_0]:",40), kapa
elseif (eos_type.eq.eos_sl) then
  write(iow,'(3X,A45)') adjl("Computation of polymer  mass density from SL EoS",45)
  write(*  ,'(3X,A45)') adjl("Computation of polymer  mass density from SL EoS",45)
  V_star             = boltz_const_Joule_K * T_star / P_star
  T_tilde            = temperature / T_star
  P_tilde            = pressure / P_star
  rsl_N              = (massOfMonomer * P_star) / (rho_star * kg_m3_to_gr_m3 * boltz_const_Joule_molK * T_star)
  rho_tilde_bulk     = eos_rho_tilde_0(T_tilde, P_tilde, rsl_N*lengthMatrix)
  massBulkDensity    = rho_tilde_bulk * rho_star
  molarBulkDensity   = massBulkDensity / massOfMonomer * gr_cm3_to_kg_m3
  segmentBulkDensity = molarBulkDensity * N_avog

  write(iow,'(3X,A45,F16.4,'' g/cm3'')')adjl('mass density was recomputed as:',45), massBulkDensity/gr_cm3_to_kg_m3
  write(*  ,'(3X,A45,F16.4,'' g/cm3'')')adjl('mass density was recomputed as:',45), massBulkDensity/gr_cm3_to_kg_m3

  sgtParam = 2.0d0 * P_star * rsl_N**2.0d0 * V_star**(8.0d0/3.0d0) * sgtParamTilde
endif
!------------------------------------------------------------------------------------------------------!
end subroutine init_scf_params
