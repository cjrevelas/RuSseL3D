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

    !calculate the radius of gyration of grafted chains
    write(iow,'(3x,A40,E16.9,A10)')adjl('Squared gyration radius of grafted ch:',40)&
                                                  &,Rg2_per_mon_gr*chainlen_gr,' Angstroms'
    write(6  ,'(3x,A40,E16.9,A10)')adjl('Squared gyration radius of grafted ch:',40)&
                                                  &,Rg2_per_mon_gr*chainlen_gr,' Angstroms'

    chainlen_matrix_max = max(chainlen_matrix, chainlen_gr)
else
    chainlen_matrix_max = chainlen_matrix
endif

ds_ave_matrix_ed   = chainlen_matrix_max / dble(ns_matrix_ed)
ds_ave_matrix_conv = chainlen_matrix / dble(ns_matrix_conv)
write(iow,'(3x,A40,E16.9)')adjl('Length of a matrix segment for Edwards:',40),ds_ave_matrix_ed
write(6  ,'(3x,A40,E16.9)')adjl('Length of a matrix segment for Edwards:',40),ds_ave_matrix_ed
write(iow,'(3x,A40,E16.9)')adjl('Length of a matrix segment for convol:',40) ,ds_ave_matrix_conv
write(6  ,'(3x,A40,E16.9)')adjl('Length of a matrix segment for convol:',40) ,ds_ave_matrix_conv

!calculate the radius of gyration of matrix chains
write(iow,'(3x,A40,E16.9,A10)')adjl('Squared gyration radius of matrix ch:',40)&
                                   &,Rg2_per_mon_matrix * chainlen_matrix,' Angstroms'
write(6  ,'(3x,A40,E16.9,A10)')adjl('Squared gyration radius of matrix ch:',40)&
                                   &,Rg2_per_mon_matrix * chainlen_matrix,' Angstroms'

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
