subroutine energies(qm_interp_mg, wa, Ufield, phia_mx, phia_gr, part_func, num_gpoints, gpid, free_energy)
!-------------------------------------------------------------------------------------------------!
use parser_vars
use geometry
use constants
use iofiles
!-------------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: num_gpoints
integer, intent(in), dimension(num_gpoints) :: gpid
integer                                     :: k1, gnode_id

real(8), intent(in), dimension(numnp,ns_gr_conv+1) :: qm_interp_mg
real(8), intent(in), dimension(numnp)              :: wa, Ufield, phia_mx, phia_gr
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

do k1 = 1, numnp
   dterm1(k1) = (1.d0 - phia_mx(k1) - phia_gr(k1))**2.d0             !blue: particle-particle interaction
   dterm2(k1) = -(wa(k1) - Ufield(k1)) * (phia_mx(k1) + phia_gr(k1)) !pink: rho-w interaction
enddo

call spat_3d(dterm1, term1, Q, vol)
call spat_3d(dterm2, term2, Q, vol)

term1 = term1 * 1.0d-30 * 0.5d0 / kappa_T
term2 = term2 * 1.0d-30 * rho_0 * boltz_const_Joule_molK * Temp
term3 = rho_0 * 1.0d-30 * vol   * boltz_const_Joule_molK * Temp * (1.d00 - part_func) / chainlen_matrix !red: entropy of matrix chains

do k1 = 1, num_gpoints
    gnode_id = gpid(k1)
    r_gpoint = min(abs(xc(3,gnode_id) - box_lo(3)), abs(box_hi(3) - xc(3,gnode_id)))

    term4_gpoint(k1) = - boltz_const_Joule_K * Temp * log(qm_interp_mg(gnode_id,ns_gr_conv+1)) !green: entropy of grafted chains
    term4            = term4 + term4_gpoint(k1) 
    
    term4_norm_gpoint(k1) = - boltz_const_Joule_K * Temp * log(r_ref / r_gpoint) 
    term4_norm            = term4_norm + term4_norm_gpoint(k1)
enddo

term1      = term1      * 1.d03 / (interf_area*1.d-20)
term2      = term2      * 1.d03 / (interf_area*1.d-20)
term3      = term3      * 1.d03 / (interf_area*1.d-20)
term4      = term4      * 1.d03 / (interf_area*1.d-20)
term4_norm = term4_norm * 1.d03 / (interf_area*1.d-20)

free_energy = term1 + term2 + term3 + term4 + term4_norm 

open(unit=837, file = energy_terms)
write(837,'(6(A19,1X))')      "term1", "term2", "term3", "term4", "term4_norm", "free_energy"
write(837,'(6(E19.9e2,1X))')   term1,   term2,   term3,   term4,   term4_norm,   free_energy

if (num_gpoints.ne.0) then
    write(837,'(A10,3(1X,A19))')  "id", "qm(ns)", "term4_gpoint", "term4_norm_gpoint"
    do k1 = 1, num_gpoints
        gnode_id   = gpid(k1)
        write(837,'(I10,3(1X,E19.9e2))') gnode_id, qm_interp_mg(gnode_id,ns_gr_conv+1), term4_gpoint(k1), term4_norm_gpoint(k1)
    enddo
endif

close(837)

return
!-------------------------------------------------------------------------------------------------!
end subroutine energies
