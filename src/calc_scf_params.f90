subroutine calc_scf_params
!------------------------------------------------------------------------------------------------------!
use parser_vars
use constants
use write_helper
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
!set the value of Hamaker constants
Aps   = Aps * 1.e-20
Asio2 = Asio2 * 1.e-20

write(iow,'(/''*Initialization of usefull quantities..'')')
write(6  ,'(/''*Initialization of usefull quantities..'')')

chainlen_free = ds_ave_free * dble(ns_free)
write(iow,'(3x,A40,E16.9)')adjl('Length of free chains:',40),chainlen_free
write(6  ,'(3x,A40,E16.9)')adjl('Length of free chains:',40),chainlen_free

!calculate the radius of gyration of free and grafted chains
write(iow,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of free ch:',40),Rg2_per_mon_free * chainlen_free,' Angstroms'
write(6  ,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of free ch:',40),Rg2_per_mon_free * chainlen_free,' Angstroms'

if (use_grafted.eq.1) then
    chainlen_gr = ds_ave_gr * dble(ns_gr)
    write(iow,'(3x,A40,E16.9)')adjl('Length of grafted chains:',40),chainlen_gr
    write(6  ,'(3x,A40,E16.9)')adjl('Length of grafted chains:',40),chainlen_gr

    write(iow,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of grafted ch:',40),Rg2_per_mon_gr*chainlen_gr,' Angstroms'
    write(6  ,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of grafted ch:',40),Rg2_per_mon_gr*chainlen_gr,' Angstroms'

    ns_free_max  = max(ns_free, ns_gr)
else
    ns_free_max = ns_free
endif

!calculate segment density in the bulk rho_0 in mol_segments/m3
rho_0 = massden/mon_mass*1.d06
write(iow,'(3x,A40,E16.9,A8)')adjl('Segment density in bulk (rho):',40),rho_0,' mol/m^3'
write(6  ,'(3x,A40,E16.9,A8)')adjl('Segment density in bulk (rho):',40),rho_0,' mol/m^3'

!calculate Helfand parameter kapa
kapa = 1.d0 /(kappa_T * boltz_const_Joule_molK * Temp * rho_0)
write(iow,'(3x,A40,E16.9)')adjl('kapa = 1/[k_T k_B T rho_0]:',40),kapa
write(6  ,'(3x,A40,E16.9)')adjl('kapa = 1/[k_T k_B T rho_0]:',40),kapa
!------------------------------------------------------------------------------------------------------!
end subroutine calc_scf_params
