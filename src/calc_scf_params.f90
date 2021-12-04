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

ds_ave_free = chainlen_free / dble(ns_free_ed)
write(iow,'(3x,A40,E16.9)')adjl('Length of a free segment:',40),ds_ave_free
write(6  ,'(3x,A40,E16.9)')adjl('Length of a free segment:',40),ds_ave_free

!calculate the radius of gyration of free and grafted chains
write(iow,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of free ch:',40),Rg2_per_mon_free * chainlen_free,' Angstroms'
write(6  ,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of free ch:',40),Rg2_per_mon_free * chainlen_free,' Angstroms'

if (use_grafted.eq.1) then
    ds_ave_gr = chainlen_gr / dble(ns_gr_ed)

    if (ds_ave_gr.ne.ds_ave_free) then
        write(ERROR_MESSAGE,'(''Lengths of grafted and free segments are not equal: '',2(E16.9))') ds_ave_gr, ds_ave_free
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    write(iow,'(3x,A40,E16.9)')adjl('Length of a grafted segment:',40),ds_ave_gr
    write(6  ,'(3x,A40,E16.9)')adjl('Length of a grafted segment:',40),ds_ave_gr

    write(iow,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of grafted ch:',40),Rg2_per_mon_gr*chainlen_gr,' Angstroms'
    write(6  ,'(3x,A40,E16.9,A10)')adjl('Square of radius of Gyration of grafted ch:',40),Rg2_per_mon_gr*chainlen_gr,' Angstroms'

    ns_free_max  = max(ns_free_ed, ns_gr_ed)
else
    ns_free_max = ns_free_ed
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
