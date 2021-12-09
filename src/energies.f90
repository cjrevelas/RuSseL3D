subroutine energies(qmx_interp_mg, phi_total, part_func, num_gpoints, gpid, free_energy)
!-------------------------------------------------------------------------------------------------!
use eos,         only: eos_ff, eos_df_drho
use parser_vars, only: ns_gr_conv, chainlen_mx, interf_area, rho_mol_bulk, temp
use geometry,    only: numnp, box_lo, box_hi, xc
use constants,   only: n_avog, boltz_const_Joule_molK, boltz_const_Joule_K
use iofiles,     only: energy_terms
!-------------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: num_gpoints
integer, intent(in), dimension(num_gpoints) :: gpid
integer                                     :: kk, gnode_id

real(8), intent(in), dimension(ns_gr_conv+1,numnp) :: qmx_interp_mg
real(8), intent(in), dimension(numnp)              :: phi_total
real(8), intent(in)                                :: part_func
real(8), intent(out)                               :: free_energy
real(8)                                            :: Q=0.d0, vol=0.d0
real(8)                                            :: term1=0.d0, term2=0.d0, term3=0.d0, term4=0.d0, term4_norm
real(8)                                            :: r_ref = 0.4d0, r_gpoint =0.d0
real(8), dimension(num_gpoints)                    :: term4_norm_gpoint, term4_gpoint
real(8), dimension(numnp)                          :: dterm1, dterm2
!-------------------------------------------------------------------------------------------------!
term4      = 0.d0
term4_norm = 0.d0
dterm1     = 0.d0
dterm2     = 0.d0

do kk = 1, numnp
   dterm1(kk) = eos_ff(phi_total(kk)) - eos_ff(1.d0)
   dterm2(kk) = -phi_total(kk) * eos_df_drho(phi_total(kk)) + 1.d0*eos_df_drho(1.d0)
enddo

call spat_3d(dterm1, term1, Q, vol)
call spat_3d(dterm2, term2, Q, vol)

term1 = term1 * 1.0d-30
term2 = term2 * 1.0d-30 * rho_mol_bulk * n_avog
term3 = rho_mol_bulk * 1.0d-30 * vol * boltz_const_Joule_molK * Temp * (1.d00 - part_func) / chainlen_mx

do kk = 1, num_gpoints
    gnode_id = gpid(kk)
    r_gpoint = MIN(ABS(xc(3,gnode_id) - box_lo(3)), ABS(box_hi(3) - xc(3,gnode_id)))

    term4_gpoint(kk) = - boltz_const_Joule_K * Temp * LOG(qmx_interp_mg(ns_gr_conv+1,gnode_id))
    term4            = term4 + term4_gpoint(kk)

    term4_norm_gpoint(kk) = - boltz_const_Joule_K * Temp * LOG(r_ref / r_gpoint)
    term4_norm            = term4_norm + term4_norm_gpoint(kk)
enddo

term1      = term1      * 1.d03 / (interf_area*1.d-20)
term2      = term2      * 1.d03 / (interf_area*1.d-20)
term3      = term3      * 1.d03 / (interf_area*1.d-20)
term4      = term4      * 1.d03 / (interf_area*1.d-20)
term4_norm = term4_norm * 1.d03 / (interf_area*1.d-20)

free_energy = term1 + term2 + term3 + term4 + term4_norm

open(unit=837, file = energy_terms)
write(837,'(A14,3A20,A22,A21)')   "term1", "term2", "term3", "term4", "term4_norm", "free_energy"
write(837,'(6(E19.9E2,1X))')   term1,   term2,   term3,   term4,   term4_norm,   free_energy

if (num_gpoints.ne.0) then
    write(837,*)
    write(837,'(A9,A17,A23,A22)')  "id", "qmx(ns)", "term4_gpoint", "term4_norm_gpoint"
    do kk = 1, num_gpoints
        gnode_id   = gpid(kk)
        write(837,'(I10,3(1X,E19.9e2))') gnode_id, qmx_interp_mg(ns_gr_conv+1,gnode_id), term4_gpoint(kk), term4_norm_gpoint(kk)
    enddo
endif

close(837)

return
!-------------------------------------------------------------------------------------------------!
end subroutine energies
