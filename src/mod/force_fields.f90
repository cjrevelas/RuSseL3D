module force_fields
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
contains

subroutine hamaker_sphere_plate(h_12, r_pol, sigma1, sigma2, A1, A2, Urep, Uatt)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: h_12, r_pol, sigma1, sigma2, A1, A2
real(8), intent(out) :: Urep, Uatt
real(8)              :: sigma = 0.d0, hamaker_constant = 0.d0
!------------------------------------------------------------------------------------------------------!
!***SI units are used everywhere***
Urep = 0.d0
Uatt = 0.d0

hamaker_constant = sqrt(A1*A2)

sigma = (sigma1 + sigma2)/2.d0

!Ruckenstein and Priere, 1976a
!Feke et al., 1984
Uatt = - (hamaker_constant/6.d0) * (1.d0/(h_12/r_pol) + 1.d0/(2.d0+h_12/r_pol) + log((h_12/r_pol)/(2.d0+h_12/r_pol)))

Urep = (hamaker_constant/7.56d3)*(sigma/r_pol)**6.d0*((8.d0 + h_12/r_pol)/(2.d0 + h_12/r_pol)**7.d0 &
                                                           & + (6.d0 - h_12/r_pol)/(h_12/r_pol)**7.d0)
return
!------------------------------------------------------------------------------------------------------!
end subroutine hamaker_sphere_plate



subroutine hamaker_sphere_sphere(h_12, r1, r2, sigma1, sigma2, A1, A2, Urep, Uatt)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: h_12, r1, r2, sigma1, sigma2, A1, A2
real(8), intent(out) :: Urep, Uatt
real(8)              :: r12 = 0.d0, sigma = 0.d0, hamaker_constant = 0.d0
!------------------------------------------------------------------------------------------------------!
!***SI units are used everywhere***
Urep = 0.d0
Uatt = 0.d0

hamaker_constant = sqrt(A1*A2)

sigma = (sigma1 + sigma2)/2.d0

r12 = r1 + r2 + h_12  !=r_centers

Uatt = -hamaker_constant/(6.d00) * (                      &
       (2*r1*r2)/(r12**2.-(r1+r2)**2.) +                  &
       (2*r1*r2)/(r12**2.-(r1-r2)**2.) +                  &
        log((r12**2.-(r1+r2)**2.)/(r12**2.-(r1-r2)**2.))  )

Urep = hamaker_constant*sigma**6./(37800.*r12)*  (                           &
       ( r12**2.-7.*r12*(r1+r2)+6*(r1**2.+7.*r1*r2+r2**2.))/(r12-r1-r2)**7.  &
      +( r12**2.+7.*r12*(r1+r2)+6*(r1**2.+7.*r1*r2+r2**2.))/(r12+r1+r2)**7.  &
      -( r12**2.+7.*r12*(r1-r2)+6*(r1**2.-7.*r1*r2+r2**2.))/(r12+r1-r2)**7.  &
      -( r12**2.-7.*r12*(r1-r2)+6*(r1**2.-7.*r1*r2+r2**2.))/(r12-r1+r2)**7.  )
return
!------------------------------------------------------------------------------------------------------!
end subroutine hamaker_sphere_sphere

!--------------------------------------------------------------------!
end module force_fields
