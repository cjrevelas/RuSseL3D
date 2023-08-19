!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine InitScfParams()
!------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: lengthMatrix, iow, sgtParam, sgtParamTilde, massDensity, massOfMonomer, &
                            massBulkDensity, molarBulkDensity, segmentBulkDensity, temperature,     &
                            pressure, matrixExist, lengthMatrix, lengthGrafted
use constants_mod,    only: boltz_const_Joule_K, boltz_const_Joule_molK, gr_cm3_to_kg_m3, N_avog, &
                            kg_m3_to_gr_m3, m3_to_cm3
use eos_mod,          only: eosType, kapa, volStar, pressStar, tempStar, rhoStar, rhoTildeBulk, &
                            pressTilde, tempTilde, rslN, helfandCompressibility, eosRhoTildeZero
use flags_mod,        only: eosHelfand, eosSanchezLacombe
use write_helper_mod, only: adjl
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8) :: lengthBulk = 0.0d0, sanchezLacombeCompressibility = 0.0d0, aux = 0.0d0
!------------------------------------------------------------------------------------------------------!
if (eosType.eq.eosHelfand) then
  molarBulkDensity   = massDensity / massOfMonomer * m3_to_cm3
  segmentBulkDensity = molarBulkDensity * N_avog

  write(iow,'(3X,A40,E16.9,A10)')adjl("Segment density in bulk (rho):",40), molarBulkDensity, " [mol/m^3]"
  write(6  ,'(3X,A40,E16.9,A10)')adjl("Segment density in bulk (rho):",40), molarBulkDensity, " [mol/m^3]"

  kapa = 1.0d0 /(helfandCompressibility * boltz_const_Joule_molK * temperature * molarBulkDensity)

  write(iow,'(3X,A40,E16.9)')adjl("kapa = 1/[k_T k_B T rho_0]:",40), kapa
  write(6  ,'(3X,A40,E16.9)')adjl("kapa = 1/[k_T k_B T rho_0]:",40), kapa
elseif (eosType.eq.eosSanchezLacombe) then
  write(iow,'(3X,A40)') adjl("Computing bulk mass density from SL EoS",40)
  write(*  ,'(3X,A40)') adjl("Computing bulk mass density from SL EoS",40)

  volStar  = boltz_const_Joule_K * tempStar / pressStar
  tempTilde = temperature / tempStar
  pressTilde = pressure / pressStar
  rslN   = (massOfMonomer * pressStar) / (rhoStar * kg_m3_to_gr_m3 * boltz_const_Joule_molK * tempStar)

  if (matrixExist.eq.1) then
    lengthBulk = lengthMatrix
  else
    lengthBulk = lengthGrafted
  endif

  write(iow,'(3X,A40,E16.9)')adjl('Chain length in the bulk:',40), lengthBulk
  write(*  ,'(3X,A40,E16.9)')adjl('Chain length in the bulk:',40), lengthBulk

  rhoTildeBulk       = eosRhoTildeZero(tempTilde, pressTilde, rslN*lengthBulk)
  massBulkDensity    = rhoTildeBulk * rhoStar
  molarBulkDensity   = massBulkDensity / massOfMonomer * gr_cm3_to_kg_m3
  segmentBulkDensity = molarBulkDensity * N_avog

  write(iow,'(3X,A40,E16.9," [g/cm3]")')adjl('Bulk mass density was recomputed as:',40), massBulkDensity/gr_cm3_to_kg_m3
  write(*  ,'(3X,A40,E16.9," [g/cm3]")')adjl('Bulk mass density was recomputed as:',40), massBulkDensity/gr_cm3_to_kg_m3

  aux = tempTilde * pressStar * rhoTildeBulk**2.0d0 * (1.0d0/(1.0d0-rhoTildeBulk) + 1.0d0/(rhoTildeBulk*rslN*lengthBulk) - 2.0d0/tempTilde)

  sanchezLacombeCompressibility = 1.0d0 / aux

  write(iow,'(3X,A40,E16.9," [Pa^-1]")')adjl('SL isothermal compressibility:',40), sanchezLacombeCompressibility
  write(*  ,'(3X,A40,E16.9," [Pa^-1]")')adjl('SL isothermal compressibility:',40), sanchezLacombeCompressibility

  sgtParam = 2.0d0 * pressStar * rslN**2.0d0 * volStar**(8.0d0/3.0d0) * sgtParamTilde
endif
!------------------------------------------------------------------------------------------------------!
end subroutine InitScfParams
