subroutine periodic_dumper(qf_final, qgr_final, phia_fr, phia_gr, wa, wa_new, wa_mix)
!-----------------------------------------------------------------------------------------------------------!
use xdata
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: k1

real(8), intent(in), dimension(numnp)      :: phia_fr, phia_gr, wa, wa_new, wa_mix
real(8), intent(in), dimension(numnp,ns+1) :: qf_final, qgr_final
!-----------------------------------------------------------------------------------------------------------!
open (unit=121, file = 'wa.out.txt', status='unknown', position='append')
open (unit=122, file = 'wa_new.out.txt', status='unknown', position='append')
open (unit=123, file = 'wa_mix.out.txt', status='unknown', position='append')
open (unit=124, file = 'rho_free.out.txt', status='unknown', position='append')
open (unit=125, file = 'rho_grafted.out.txt', status='unknown', position='append')

do k1 = 1, numnp, print_ev
   write(121,'(E19.9e3)',advance='no') wa(k1)
   write(122,'(E19.9e3)',advance='no') wa_new(k1)
   write(123,'(E19.9e3)',advance='no') wa_mix(k1)
   write(124,'(E19.9e3)',advance='no') phia_fr(k1)
   write(125,'(E19.9e3)',advance='no') phia_gr(k1)
enddo

write(121,*)
write(122,*)
write(123,*)
write(124,*)
write(125,*)

close(121)
close(122)
close(123)
close(124)
close(125)

open (unit=120, file = 'rho_reduced.out.txt')
write(120,'(A13,8(A19))') 'np','x','y','z','phi_f','phi_g','wa','wa_new','wa_mix'
do k1 = 1, numnp
   write(120,'(I13,8(E19.9e3))') k1, xc(1,k1), xc(2,k1), xc(3,k1), phia_fr(k1), &
   &                             phia_gr(k1), wa(k1), wa_new(k1), wa_mix(k1)
enddo
close(120)

!print the restricted partition functions
call qprint(qf_final,"free")

if (use_grafted.eq.1) then
   call qprint(qgr_final,"graf")
endif
!-----------------------------------------------------------------------------------------------------------!
end subroutine periodic_dumper
