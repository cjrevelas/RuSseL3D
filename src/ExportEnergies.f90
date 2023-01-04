!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportEnergies(qmx_interp_mg, qgr_interp, phi_total, ww, Ufield, part_func, &
                           targetNumGraftedChains, graftPointId, free_energy)
!-------------------------------------------------------------------------------------------------!
use eos_mod,         only: eos_ff, eos_df_drho
use parser_vars_mod, only: numConvolPointsGrafted, lengthMatrix, molarBulkDensity, temperature, &
                           beta, graftPointDistance, matrixExist, graftedExist
use geometry_mod,    only: numNodes, interfaceArea
use constants_mod,   only: n_avog, boltz_const_Joule_molK, N_to_mN, A2_to_m2, A3_to_m3
use iofiles_mod,     only: IO_energies
!-------------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------------!
integer, intent(in)                                    :: targetNumGraftedChains
integer, intent(in), dimension(targetNumGraftedChains) :: graftPointId
integer                                                :: kk, gnode_id


real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qmx_interp_mg, qgr_interp
real(8), intent(in), dimension(numNodes)                          :: phi_total, ww, Ufield
real(8), intent(in)                                               :: part_func
real(8), intent(out)                                              :: free_energy
real(8)                                                           :: ComputeStretchingEnergy
real(8)                                                           :: Q=0.0d0, vol=0.0d0
real(8)                                                           :: E_eos_f=0.0d0, E_eos_dfdrho=0.0d0
real(8)                                                           :: E_field=0.0d0, E_solid=0.0d0
real(8)                                                           :: E_entropy_mx=0.0d0, E_entropy_gr=0.0d0
real(8)                                                           :: E_entropy_gr_normlz=0.0d0, E_stretching = 0.0d0
real(8)                                                           :: r_ref=0.4d0
real(8)                                                           :: wwFieldBulk=0.0d0
real(8), dimension(targetNumGraftedChains)                        :: term4_norm_gpoint, term4_gpoint
real(8), dimension(numNodes)                                      :: prof_eos_f, prof_eos_dfdrho
real(8), dimension(numNodes)                                      :: prof_field, prof_solid
!-------------------------------------------------------------------------------------------------!
E_entropy_gr        = 0.0d0
E_entropy_gr_normlz = 0.0d0
E_stretching        = 0.0d0

prof_eos_f      = 0.0d0
prof_eos_dfdrho = 0.0d0
prof_field      = 0.0d0
prof_solid      = 0.0d0

wwFieldBulk = eos_df_drho(1.0d0) * beta

do kk = 1, numNodes
  prof_eos_f(kk)      =  eos_ff(phi_total(kk))
  prof_eos_dfdrho(kk) = -phi_total(kk)*eos_df_drho(phi_total(kk))
  prof_field(kk)      = -phi_total(kk)*(ww(kk) + wwFieldBulk)
  prof_solid(kk)      =  phi_total(kk)*Ufield(kk)
enddo

if (matrixExist.eq.1) then
  do kk = 1, numNodes
    prof_eos_f(kk)      = prof_eos_f(kk) - eos_ff(1.0d0)
    prof_eos_dfdrho(kk) = prof_eos_dfdrho(kk) + 1.0d0*eos_df_drho(1.0d0)
    prof_field(kk)      = prof_field(kk) - (-1.0d0 * wwFieldBulk)
  enddo
endif

call FemIntegration(prof_eos_f, E_eos_f, Q, vol)
call FemIntegration(prof_eos_dfdrho, E_eos_dfdrho, Q, vol)
call FemIntegration(prof_field, E_field, Q, vol)
call FemIntegration(prof_solid, E_solid, Q, vol)

E_eos_f      = E_eos_f      * A3_to_m3
E_eos_dfdrho = E_eos_dfdrho * A3_to_m3 * molarBulkDensity * n_avog
E_field      = E_field      * A3_to_m3 * molarBulkDensity * n_avog / beta
E_solid      = E_solid      * A3_to_m3 * molarBulkDensity * n_avog / beta

if (matrixExist.eq.1) E_entropy_mx = molarBulkDensity * A3_to_m3 * vol * boltz_const_Joule_molK * temperature * (1.0d0 - part_func) / lengthMatrix

if (graftedExist.eq.1) then
  do kk = 1, targetNumGraftedChains
    gnode_id              =  graftPointId(kk)
    term4_gpoint(kk)      = -LOG(qmx_interp_mg(numConvolPointsGrafted+1,gnode_id)) / beta
    term4_norm_gpoint(kk) = -LOG(r_ref/graftPointDistance) / beta  ! TODO: this needs revision in the case of multiple solid surfaces
    E_entropy_gr          =  E_entropy_gr        + term4_gpoint(kk)
    E_entropy_gr_normlz   =  E_entropy_gr_normlz + term4_norm_gpoint(kk)
    E_stretching          =  E_stretching        + ComputeStretchingEnergy(gnode_id, qmx_interp_mg, qgr_interp)
  enddo
endif

E_eos_f             = E_eos_f             * N_to_mN / (interfaceArea()*A2_to_m2)
E_eos_dfdrho        = E_eos_dfdrho        * N_to_mN / (interfaceArea()*A2_to_m2)
E_field             = E_field             * N_to_mN / (interfaceArea()*A2_to_m2)
E_solid             = E_solid             * N_to_mN / (interfaceArea()*A2_to_m2)

if (matrixExist.eq.1) E_entropy_mx        = E_entropy_mx        * N_to_mN / (interfaceArea()*A2_to_m2)

if (graftedExist.eq.1) then
  E_entropy_gr        = E_entropy_gr        * N_to_mN / (interfaceArea()*A2_to_m2)
  E_entropy_gr_normlz = E_entropy_gr_normlz * N_to_mN / (interfaceArea()*A2_to_m2)
  E_stretching        = E_stretching        * N_to_mN / (interfaceArea()*A2_to_m2)
endif

free_energy = E_eos_f + E_eos_dfdrho + E_entropy_mx + E_entropy_gr + E_entropy_gr_normlz

open(unit=837, file = IO_energies)
write(837,'(9(2X,A16))')     "E_ff", "E_dfdrho", "E_entropy_mx", "E_entr_gr", "E_entr_gr_norm", "E_total", "E_stretch", "E_w", "E_Us"
write(837,'(9(2X,E16.9E2))') E_eos_f, E_eos_dfdrho, E_entropy_mx, E_entropy_gr, E_entropy_gr_normlz, free_energy, E_stretching, E_field, E_solid

if (graftedExist.eq.1) then
  if (targetNumGraftedChains.ne.0) then
    write(837,*)
    write(837,'(4(2X,A16))')  "id", "qmx(ns)", "E_entr_gp", "E_entr_gp_norm"
    do kk = 1, targetNumGraftedChains
      gnode_id   = graftPointId(kk)
      write(837,'(2X,I16,3(2X,E16.9E2))') gnode_id, qmx_interp_mg(numConvolPointsGrafted+1,gnode_id), term4_gpoint(kk), term4_norm_gpoint(kk)
    enddo
  endif
endif

close(837)

return
!-------------------------------------------------------------------------------------------------!
end subroutine ExportEnergies
