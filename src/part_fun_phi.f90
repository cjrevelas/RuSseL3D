subroutine part_fun_phi(q_final, phia, part_func, nch_per_area, coef)    !CJR
!-------------------------------------------------------------------------------------------! 
use xdata
!-------------------------------------------------------------------------------------------!     
implicit none
!-------------------------------------------------------------------------------------------!    
integer :: k1                                            !CJR

real(8), intent(in), dimension(numnp,ns+1) :: q_final    !CJR
real(8), intent(in), dimension(numnp)      :: phia       !CJR
real(8), intent(out) :: part_func, nch_per_area, coef    !CJR
real(8) :: sum_f, Q
real(8), dimension(numnp) :: q_last
!-------------------------------------------------------------------------------------------! 
sum_f = 0.d0

do k1 = 1, numnp
    q_last(k1) = q_final(k1,ns+1)      !CJR
enddo

call spat_3d(q_last, sum_f, Q)

part_func = Q



sum_f = 0.d0

call spat_3d(phia, sum_f, Q)

nch_per_area = sum_f*1.0d-30*rho_0/chainlen

coef = part_func/nch_per_area/(chainlen/(rho_0*volume*1.d-30))

return
!-------------------------------------------------------------------------------------------!
end subroutine part_fun_phi
