subroutine calc_scf_params
!------------------------------------------------------------------------------------------------------!
use parser_vars
use constants
use write_helper
use error_handing
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
write(iow,'(/''*Initialization of usefull quantities..'')')
write(6  ,'(/''*Initialization of usefull quantities..'')')

if (use_grafted.eq.1) then
    ds_ave_gr_ed   = chainlen_gr / dble(ns_gr_ed)
    ds_ave_gr_conv = chainlen_gr / dble(ns_gr_conv)

    write(iow,'(3x,A40,E16.9)')adjl('Length of a grafted segment for Edwards:',40),ds_ave_gr_ed
    write(6  ,'(3x,A40,E16.9)')adjl('Length of a grafted segment for Edwards:',40),ds_ave_gr_ed
    write(iow,'(3x,A40,E16.9)')adjl('Length of a grafted segment for convol:',40) ,ds_ave_gr_conv
    write(6  ,'(3x,A40,E16.9)')adjl('Length of a grafted segment for convol:',40) ,ds_ave_gr_conv

    write(iow,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of grafted ch:',40),Rg2_per_mon_gr*chainlen_gr,' Angstroms'
    write(6  ,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of grafted ch:',40),Rg2_per_mon_gr*chainlen_gr,' Angstroms'

    chainlen_free_max = max(chainlen_free, chainlen_gr)
else
    chainlen_free_max = chainlen_free
endif

ds_ave_free_ed   = chainlen_free_max / dble(ns_free_ed)
ds_ave_free_conv = chainlen_free / dble(ns_free_conv)
write(iow,'(3x,A40,E16.9)')adjl('Length of a free segment for Edwards:',40),ds_ave_free_ed
write(6  ,'(3x,A40,E16.9)')adjl('Length of a free segment for Edwards:',40),ds_ave_free_ed
write(iow,'(3x,A40,E16.9)')adjl('Length of a free segment for convol:',40) ,ds_ave_free_conv
write(6  ,'(3x,A40,E16.9)')adjl('Length of a free segment for convol:',40) ,ds_ave_free_conv

!calculate the radius of gyration of free and grafted chains
write(iow,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of free ch:',40),Rg2_per_mon_free * chainlen_free,' Angstroms'
write(6  ,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of free ch:',40),Rg2_per_mon_free * chainlen_free,' Angstroms'

!calculate segment density in the bulk, rho_0, in mol_segments/m3
rho_0 = massden/mon_mass*1.d06
write(iow,'(3x,A40,E16.9,A8)')adjl('Segment density in bulk (rho):',40),rho_0,' mol/m^3'
write(6  ,'(3x,A40,E16.9,A8)')adjl('Segment density in bulk (rho):',40),rho_0,' mol/m^3'

!calculate Helfand parameter kapa
kapa = 1.d0 /(kappa_T * boltz_const_Joule_molK * Temp * rho_0)
write(iow,'(3x,A40,E16.9)')adjl('kapa = 1/[k_T k_B T rho_0]:',40),kapa
write(6  ,'(3x,A40,E16.9)')adjl('kapa = 1/[k_T k_B T rho_0]:',40),kapa

!------------------------------------------------------------------------------------------------------!
end subroutine calc_scf_params
