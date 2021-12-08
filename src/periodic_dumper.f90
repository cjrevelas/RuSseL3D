subroutine periodic_dumper(qm_final, qgr_final, qm_interp_mm, qm_interp_mg, qgr_interp, phia_mx, phia_gr, wa, wa_new, wa_mix)
!-----------------------------------------------------------------------------------------------------------!
use parser_vars
use geometry
use iofiles
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: k1

real(8), intent(in), dimension(numnp)                  :: phia_mx, phia_gr, wa, wa_new, wa_mix
real(8), intent(in), dimension(ns_matrix_ed+1,numnp)   :: qm_final
real(8), intent(in), dimension(ns_matrix_conv+1,numnp) :: qm_interp_mm
real(8), intent(in), dimension(ns_gr_conv+1,numnp)     :: qm_interp_mg
real(8), intent(in), dimension(ns_gr_ed+1,numnp)       :: qgr_final
real(8), intent(in), dimension(ns_gr_conv+1,numnp)     :: qgr_interp
!-----------------------------------------------------------------------------------------------------------!
open (unit=120, file = profiles)
write(120,'(A13,A11,2A19,2A21,A15,A22,A19)') "np", "x", "y", "z", "phi_m", "phi_g", "wa", "wa_new", "wa_mix"
do k1 = 1, numnp
   write(120,'(I13,8(E19.9E3))') k1, xc(1,k1), xc(2,k1), xc(3,k1), phia_mx(k1), &
   &                             phia_gr(k1), wa(k1), wa_new(k1), wa_mix(k1)
enddo
close(120)

!print the restricted partition functions
call qprint(ns_matrix_ed, qm_final, "mtrx")
#ifdef DEBUG_OUTPUTS
call qprint(ns_matrix_conv, qm_interp_mm, "mxcf")
call qprint(ns_gr_conv, qm_interp_mg, "mxcg")
#endif

if (grafted_exist.eq.1) then
   call qprint(ns_gr_ed, qgr_final, "graf")
#ifdef DEBUG_OUTPUTS
   call qprint(ns_gr_conv, qgr_interp, "grcv")
#endif
endif
!-----------------------------------------------------------------------------------------------------------!
end subroutine periodic_dumper
