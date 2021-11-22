  subroutine part_fun_phi
   
    use xdata
    use mdata
    
    implicit none
    
    real(8):: sum_f,q_last(numnp),phi_int(numnp),Q
                   

!===================================================================!
!                                                                   !
!           Partition fun!tion Calculation- Save Field              !
!                                                                   !
!===================================================================!
!
!                          lx
!                           /
!                          |
!                      1   |
!                 Q = ---  | q(x,1) dx
!                      lx  |
!                          |
!                          /
!                        0
!
!-------------------------------------------------------------------!

        

         sum_f = 0.d0
         do k1 =1,numnp
              q_last(k1)=qf_final(k1,ns)
         end do
         
         call spat_3d(q_last,sum_f,Q)
		 part_func=Q
   

          
 

          
!===================================================================!
!
!           Reduced Density (phi) 
!
!===================================================================!

!
!              1
!             /
!            |
!            |
!  phi(x) =  | q(x,s)  q(x,1-s) ds
!            |
!            |
!           /
!          0
!
!........Estimate the reduced densities phi(x) = rho(x)/rho_0
!........Estimate the free energy,
!........using the Simpsons formula
!
!
!........ phi(x) = Int_{from 0 to 1} q(x,s) q(x,1-s) ds
!
!
!------------------------------------------------------------
!
!
!............Here we will use the Simpson formula
!............for integration over s
!
!
!
         do k1 = 1,numnp
              sum = 0.d0
              do time_step = 1, ns+1
                   sum = sum +	koeff(time_step)*qf_final(k1,time_step)* qf_final(k1,ns+2-time_step)*ds
              end do
              phia_new(k1) = sum
         end do

         sum_f = 0.d0
         call spat_3d(phia_new,sum_f,Q)
  
         nch_per_area = sum_f * 1.0d-10 * rho_0 / dfloat(chainlen)
         coef= part_func /nch_per_area/(dfloat(chainlen) / (rho_0 * lx  * 1.d-10))
! nch_per_area = sum * 1.0d-30 * N_avog*rho_0 / dfloat(chainlen)

    return 
end
