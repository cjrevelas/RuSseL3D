subroutine init_files
!------------------------------------------------------------------------------------------------------!
use parser_vars
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
if (init_iter.eq.0) then
    open(unit=121, file = 'wa.out.txt', status='replace')
    close(121)
    open(unit=121, file = 'wa_new.out.txt', status='replace')
    close(121)
    open(unit=121, file = 'wa_mix.out.txt', status='replace')
    close(121)
    open(unit=121, file = 'rho.out.txt', status='replace')
    close(121)
    open(unit=121, file = 'energy_terms.out.txt', status='replace')
    close(121)
endif
!------------------------------------------------------------------------------------------------------!
end subroutine init_files
