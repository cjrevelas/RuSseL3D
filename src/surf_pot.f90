function surf_pot(h12)
!--------------------------------------------------------------------!     
use xdata
use constants
!--------------------------------------------------------------------! 
implicit none
!--------------------------------------------------------------------! 
integer :: simple

real(8) :: surf_pot
real(8) :: sigma

real(8) :: h12, Urep, Uattr

real(8) :: a1, a2, r12, h_12

real(8) :: hamaker_constant

real(8) :: number_density
!--------------------------------------------------------------------!     
simple=0

open(unit=43, file = 'Hamaker_potential.out.txt')
 
hamaker_constant = sqrt(Aps*Asio2)

sigma = (sigma1 + sigma2)/2. * 1.e-10

h_12 = h12/1.e+10

number_density = rho_0*avogadro_constant

a1 = (3./4./pi/number_density)**(1./3.)

a2 = sphere_radius/1.e+10

r12 = a1 + a2 + h_12

surf_pot = 0.d0

if (simple.eq.0)then   
    if (h12.gt.1.e-30) then
           
        Uattr = - hamaker_constant/(6.d00) * (         &

                  (2*a1*a2)/(r12**2.-(a1+a2)**2.) +    &

                  (2*a1*a2)/(r12**2.-(a1-a2)**2.) +    &

                  log((r12**2.-(a1+a2)**2.)/(r12**2.-(a1-a2)**2.))  ) 


        Urep = hamaker_constant*sigma**6./(37800.*r12)*                                 &
                  (                                                                     &
                  ( r12**2.-7.*r12*(a1+a2)+6*(a1**2.+7.*a1*a2+a2**2.))/(r12-a1-a2)**7.  &
        
                 +( r12**2.+7.*r12*(a1+a2)+6*(a1**2.+7.*a1*a2+a2**2.))/(r12+a1+a2)**7.  &
           
                 -( r12**2.+7.*r12*(a1-a2)+6*(a1**2.-7.*a1*a2+a2**2.))/(r12+a1-a2)**7.  &
            
                 -( r12**2.-7.*r12*(a1-a2)+6*(a1**2.-7.*a1*a2+a2**2.))/(r12-a1+a2)**7.  )    

        surf_pot = Uattr + Urep
    endif   
else  
    if (h12.gt.1.e-30) then

        Uattr = - hamaker_constant/(6.d00*pi*h_12*h_12*h_12*number_density)
  
        Urep  = hamaker_constant*sigma**6. &

                /(45.d00*pi*number_density*h_12**9.)

        surf_pot = Uattr + Urep
    endif  
endif 

surf_pot = surf_pot/(boltz_const_Joule_K*Temp)
write(43,'(E16.9,2X,E16.9)') h12, surf_pot

surf_pot = surf_pot * chainlen

if (surf_pot.gt.0.) then
    surf_pot = 0.d0
endif

return 
!--------------------------------------------------------------------! 
end function surf_pot
   

