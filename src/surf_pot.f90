subroutine surf_pot(Temp, h12, sphere_radius, rho_0, sigma1, sigma2, Aps, Asio2, r12, Urep, Uatt, Utot)
!------------------------------------------------------------------------------------------------------!
use error_handing
use constants
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: Temp, h12, sphere_radius, rho_0, sigma1, sigma2, Aps, Asio2
real(8), intent(out) :: r12, Urep, Uatt, Utot
real(8)              :: sigma, hamaker_constant, number_density
real(8)              :: a1, a2, h_12
!------------------------------------------------------------------------------------------------------!
!***SI units are used everywhere***
r12  = 0.d0
Urep = 0.d0
Uatt = 0.d0
Utot = 0.d0

hamaker_constant = sqrt(Aps*Asio2)

sigma = (sigma1 + sigma2)/2. * 1.e-10

h_12 = h12/1.e+10

number_density = rho_0*avogadro_constant

a1  = (3./4./pi/number_density)**(1./3.)
a2  = sphere_radius/1.e+10
r12 = a1 + a2 + h_12

if (h12.gt.1.e-10) then
    Uatt = -hamaker_constant/(6.d00) * (                            &
           (2*a1*a2)/(r12**2.-(a1+a2)**2.) +                    &
           (2*a1*a2)/(r12**2.-(a1-a2)**2.) +                    &
            log((r12**2.-(a1+a2)**2.)/(r12**2.-(a1-a2)**2.))  )

    Urep = hamaker_constant*sigma**6./(37800.*r12)*  (                               &
           ( r12**2.-7.*r12*(a1+a2)+6*(a1**2.+7.*a1*a2+a2**2.))/(r12-a1-a2)**7.  &
          +( r12**2.+7.*r12*(a1+a2)+6*(a1**2.+7.*a1*a2+a2**2.))/(r12+a1+a2)**7.  &
          -( r12**2.+7.*r12*(a1-a2)+6*(a1**2.-7.*a1*a2+a2**2.))/(r12+a1-a2)**7.  &
          -( r12**2.-7.*r12*(a1-a2)+6*(a1**2.-7.*a1*a2+a2**2.))/(r12-a1+a2)**7.  )

    Utot = Uatt + Urep
endif

Uatt = Uatt/(boltz_const_Joule_K*Temp)
Urep = Urep/(boltz_const_Joule_K*Temp)
Utot = Utot/(boltz_const_Joule_K*Temp)

if (Utot.gt.0.) then
    Utot = 0.d0
endif

r12 = r12 * 1.e+10

return
!------------------------------------------------------------------------------------------------------!
end subroutine surf_pot
