       
function surf_pot(h12)
    
     use xdata

    implicit none

    real*8 surf_pot, U_solid
    real*8 As,eps,sigma,coeff!,layerD
    real*8 h12,Urep,Uattr
    real*8 a1,a2,r12,h_12
    real*8 hamaker_constant
    real*8 number_density!,rho_0,Temp
    integer simple!chainlen,

      
         simple=0
         open(file='Hamaker_potential.txt', unit=43)
      
     
         hamaker_constant=sqrt(Aps*Asio2)
!    :  -Aps
!      -Aps
        
         sigma=(sigma1+sigma2)/2.*1.e-10
         h_12=h12/1.e+10
         a2=sphere_radius/1.e+10
         number_density=rho_0*N_avog
         a1= (3./4./pi/number_density)**(1./3.)
         r12=a1+a2+h_12
         surf_pot=0.
       if (simple.eq.0)then
        
               if(h12.gt.1.e-30)then
           
                   Uattr=-  hamaker_constant/(6.d00)*(                &
                  (2*a1*a2)/(r12**2.-(a1+a2)**2.)+                    &
	              (2*a1*a2)/(r12**2.-(a1-a2)**2.)+                    &
	              log((r12**2.-(a1+a2)**2.)/(r12**2.-(a1-a2)**2.))      &
                   )
           
           
                    Urep=hamaker_constant*sigma**6./(37800.*r12)*                             &
                  (                                                                           &
                  ( r12**2.-7.*r12*(a1+a2)+6*(a1**2.+7.*a1*a2+a2**2.))/(r12-a1-a2)**7.        &
        
                 +( r12**2.+7.*r12*(a1+a2)+6*(a1**2.+7.*a1*a2+a2**2.))/(r12+a1+a2)**7.        &
           
                 -( r12**2.+7.*r12*(a1-a2)+6*(a1**2.-7.*a1*a2+a2**2.))/(r12+a1-a2)**7.        &
            
                 -( r12**2.-7.*r12*(a1-a2)+6*(a1**2.-7.*a1*a2+a2**2.))/(r12-a1+a2)**7.        &
                  )
                     surf_pot=Uattr+Urep
            end if 
           
           
         else 
              if(h12.gt.1.e-30)then
                   Uattr=-hamaker_constant/(6.d00*pi*h_12*h_12*h_12*number_density)
                   Urep=hamaker_constant*sigma**6.&
                        /(45.d00*pi*number_density*h_12**9.)
                   surf_pot=Uattr+Urep
              end if  
      
         end if 
           
           
           
           
       surf_pot=surf_pot/(K_boltzmann*Temp)
      write(43,'(e16.9,2x,e16.9)')h12, surf_pot
         

      surf_pot=surf_pot *dfloat (chainlen) 
       
            !!!!!!!!!!!!!!! 
 
       If(surf_pot.gt.0.)then
                      surf_pot=0.
                      endif


      






!  if (distance.le.wwidth) then
	   
!	     = - pobyV_n*Vored
	    
!     else

!         surf_pot = 0.d0

!	end if
!
	return 
!
	end
	    

