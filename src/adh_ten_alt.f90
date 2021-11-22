  subroutine adhesion_tension_alternative
    
   
    use xdata
    use mdata
   
   
    implicit none

    real*8 phi_adh(numnp),sum_f,Q
!---------------------------------------------------------------------
!**********************************************************************
!********************Grand potential calculations**********************
!**********************************************************************
!**********************************************************************
!---------------------------------------------------------------------

         sum_f = 0.d00 
		 do k1 = 1,numnp
              phi_adh(k1) = 0.5d0 * kapa * (1.d0-(phia_new(k1))**2.)
         end do
         call spat_3d(phi_adh,sum_f,Q)
!........sum_final is part of Omega N/(rho_0 k_B T a), measured in meters.
		 sum_f= sum_f * 1.d-10
         term3=xc(1,numnp)*1.0d-10 *(1.d00 - part_func)
!DEBUG
!   term4=-dble(chainlen)/rho_0*log(qf_final(1,ns))
!DEBUG
          sum_f = sum_f +term3

!........Adhesion tension, in mN/m
!........adh_ten = (gamma_s - gamma_sl)
          adh_ten_alt = - sum_f* rho_0 * k_B* Temp / dfloat(chainlen) * 1.0D+03

    return 
end
