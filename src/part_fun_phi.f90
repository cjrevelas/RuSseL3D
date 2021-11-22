subroutine part_fun_phi
!-------------------------------------------------------------------------------------------! 
use xdata
use mdata
!-------------------------------------------------------------------------------------------!     
implicit none
!-------------------------------------------------------------------------------------------!    
real(8) :: sum_f, q_last(numnp), phi_int(numnp), Q
!-------------------------------------------------------------------------------------------! 
sum_f = 0.d0

do k1 = 1, numnp
    q_last(k1) = qf_final(k1,ns)
enddo

call spat_3d(q_last, sum_f, Q)

part_func = Q

do k1 = 1, numnp
    sum = 0.d0

    do time_step = 1, ns+1
       sum = sum + koeff(time_step)*qf_final(k1,time_step)*qf_final(k1,ns+2-time_step)*ds
    enddo

    phia_new(k1) = sum
enddo

sum_f = 0.d0
call spat_3d(phia_new, sum_f, Q)
  
nch_per_area = sum_f*1.0d-10*rho_0/dfloat(chainlen)

coef= part_func/nch_per_area/(dfloat(chainlen)/(rho_0*lx*1.d-10))

return 
!-------------------------------------------------------------------------------------------!
end subroutine part_fun_phi
