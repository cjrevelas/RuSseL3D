!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_energies(qmx_interp_mg, qgr_interp, phi_total, wa, Ufield, part_func, num_gpoints, gpid, free_energy)
!-------------------------------------------------------------------------------------------------!
use eos,         only: eos_ff, eos_df_drho
use parser_vars, only: ns_gr_conv, chainlen_mx, rho_mol_bulk, temp, beta, r_gpoint
use geometry,    only: numnp, interf_area
use constants,   only: n_avog, boltz_const_Joule_molK, N_to_mN, A2_to_m2, A3_to_m3
use iofiles,     only: energy_terms
!-------------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: num_gpoints
integer, intent(in), dimension(num_gpoints) :: gpid
integer                                     :: kk, gnode_id

real(8), intent(in), dimension(ns_gr_conv+1,numnp) :: qmx_interp_mg, qgr_interp
real(8), intent(in), dimension(numnp)              :: phi_total, wa, Ufield
real(8), intent(in)                                :: part_func
real(8), intent(out)                               :: free_energy
real(8)                                            :: compute_stretching_energy
real(8)                                            :: Q=0.d0, vol=0.d0
real(8)                                            :: E_eos_f=0.d0, E_eos_dfdrho=0.d0
real(8)                                            :: E_field=0.d0, E_solid=0.d0
real(8)                                            :: E_entropy_mx=0.d0, E_entropy_gr=0.d0
real(8)                                            :: E_entropy_gr_normlz=0.d0, E_stretching = 0.d0
real(8)                                            :: r_ref=0.4d0
real(8), dimension(num_gpoints)                    :: term4_norm_gpoint, term4_gpoint
real(8), dimension(numnp)                          :: prof_eos_f, prof_eos_dfdrho
real(8), dimension(numnp)                          :: prof_field, prof_solid
!-------------------------------------------------------------------------------------------------!
E_entropy_gr        = 0.d0
E_entropy_gr_normlz = 0.d0
E_stretching        = 0.d0

prof_eos_f      = 0.d0
prof_eos_dfdrho = 0.d0
prof_field      = 0.d0
prof_solid      = 0.d0

do kk = 1, numnp
   prof_eos_f(kk)      =  eos_ff(phi_total(kk)) - eos_ff(1.d0)
   prof_eos_dfdrho(kk) = -phi_total(kk)*eos_df_drho(phi_total(kk)) + 1.d0*eos_df_drho(1.d0)
   prof_field(kk)      = -phi_total(kk)*wa(kk)
   prof_solid(kk)      =  phi_total(kk)*Ufield(kk)
enddo

call fem_integration(prof_eos_f, E_eos_f, Q, vol)
call fem_integration(prof_eos_dfdrho, E_eos_dfdrho, Q, vol)
call fem_integration(prof_field, E_field, Q, vol)
call fem_integration(prof_solid, E_solid, Q, vol)

E_eos_f      = E_eos_f      * A3_to_m3
E_eos_dfdrho = E_eos_dfdrho * A3_to_m3 * rho_mol_bulk * n_avog
E_field      = E_field      * A3_to_m3 * rho_mol_bulk * n_avog / beta
E_solid      = E_solid      * A3_to_m3 * rho_mol_bulk * n_avog / beta
E_entropy_mx = rho_mol_bulk * A3_to_m3 * vol * boltz_const_Joule_molK * Temp * (1.d00 - part_func) / chainlen_mx

do kk = 1, num_gpoints
    gnode_id              =  gpid(kk)
    term4_gpoint(kk)      = -LOG(qmx_interp_mg(ns_gr_conv+1,gnode_id)) / beta
    term4_norm_gpoint(kk) = -LOG(r_ref/r_gpoint) / beta
    E_entropy_gr          =  E_entropy_gr        + term4_gpoint(kk)
    E_entropy_gr_normlz   =  E_entropy_gr_normlz + term4_norm_gpoint(kk)
    E_stretching          =  E_stretching        + compute_stretching_energy(gnode_id, qmx_interp_mg, qgr_interp)
enddo

E_eos_f             = E_eos_f             * N_to_mN / (interf_area*A2_to_m2)
E_eos_dfdrho        = E_eos_dfdrho        * N_to_mN / (interf_area*A2_to_m2)
E_field             = E_field             * N_to_mN / (interf_area*A2_to_m2)
E_solid             = E_solid             * N_to_mN / (interf_area*A2_to_m2)
E_entropy_mx        = E_entropy_mx        * N_to_mN / (interf_area*A2_to_m2)
E_entropy_gr        = E_entropy_gr        * N_to_mN / (interf_area*A2_to_m2)
E_entropy_gr_normlz = E_entropy_gr_normlz * N_to_mN / (interf_area*A2_to_m2)
E_stretching        = E_stretching        * N_to_mN / (interf_area*A2_to_m2)

free_energy = E_eos_f + E_eos_dfdrho + E_entropy_mx + E_entropy_gr + E_entropy_gr_normlz

open(unit=837, file = energy_terms)
write(837,'(9(2X,A16))')     "E_ff", "E_dfdrho", "E_entropy_mx", "E_entr_gr", "E_entr_gr_norm", "E_total", "E_stretch", "E_w", "E_Us"
write(837,'(9(2X,E16.9E2))')E_eos_f, E_eos_dfdrho, E_entropy_mx, E_entropy_gr, E_entropy_gr_normlz, free_energy, E_stretching, E_field, E_solid

if (num_gpoints.ne.0) then
    write(837,*)
    write(837,'(4(2X,A16))')  "id", "qmx(ns)", "E_entr_gp", "E_entr_gp_norm"
    do kk = 1, num_gpoints
        gnode_id   = gpid(kk)
        write(837,'(2X,I16,3(2X,E16.9E2))') gnode_id, qmx_interp_mg(ns_gr_conv+1,gnode_id), term4_gpoint(kk), term4_norm_gpoint(kk)
    enddo
endif

close(837)

return
!-------------------------------------------------------------------------------------------------!
end subroutine compute_energies
