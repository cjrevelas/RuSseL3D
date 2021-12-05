subroutine periodic_dumper(qm_final, qgr_final, qm_interp_mm, qm_interp_mg, qgr_interp, phia_mx, phia_gr, wa, wa_new, wa_mix)
!-----------------------------------------------------------------------------------------------------------!
use parser_vars
use geometry
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: k1

real(8), intent(in), dimension(numnp)                :: phia_mx, phia_gr, wa, wa_new, wa_mix
real(8), intent(in), dimension(numnp,ns_matrix_ed+1)   :: qm_final
real(8), intent(in), dimension(numnp,ns_matrix_conv+1) :: qm_interp_mm
real(8), intent(in), dimension(numnp,ns_gr_conv+1)   :: qm_interp_mg
real(8), intent(in), dimension(numnp,ns_gr_ed+1)     :: qgr_final
real(8), intent(in), dimension(numnp,ns_gr_conv+1)   :: qgr_interp
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
   write(124,'(E19.9e3)',advance='no') phia_mx(k1)
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
   write(120,'(I13,8(E19.9e3))') k1, xc(1,k1), xc(2,k1), xc(3,k1), phia_mx(k1), &
   &                             phia_gr(k1), wa(k1), wa_new(k1), wa_mix(k1)
enddo
close(120)

!print the restricted partition functions
call qprint(ns_matrix_ed, qm_final, "free")
call qprint(ns_matrix_conv, qm_interp_mm, "frcf")
call qprint(ns_gr_conv, qm_interp_mg, "frcg")

if (use_grafted.eq.1) then
   call qprint(ns_gr_ed, qgr_final, "graf")
   call qprint(ns_gr_conv, qgr_interp, "grcv")
endif
!-----------------------------------------------------------------------------------------------------------!
end subroutine periodic_dumper
