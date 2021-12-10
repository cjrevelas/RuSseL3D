!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

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

hamaker_constant = SQRT(A1*A2)

sigma = (sigma1 + sigma2)/2.d0

! Ruckenstein and Prieve, 1976a; Feke et al., 1984
Uatt = - (hamaker_constant/6.d0) * (1.d0/(h_12/r_pol) + 1.d0/(2.d0+h_12/r_pol) + LOG((h_12/r_pol)/(2.d0+h_12/r_pol)))

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

hamaker_constant = SQRT(A1*A2)

sigma = (sigma1 + sigma2)/2.d0

r12 = r1 + r2 + h_12  !=r_centers

Uatt = -hamaker_constant/(6.d00) * (                      &
       (2*r1*r2)/(r12**2.-(r1+r2)**2.) +                  &
       (2*r1*r2)/(r12**2.-(r1-r2)**2.) +                  &
        LOG((r12**2.-(r1+r2)**2.)/(r12**2.-(r1-r2)**2.))  )

Urep = hamaker_constant*sigma**6./(37800.*r12)*  (                         &
       (r12**2.-7.*r12*(r1+r2)+6*(r1**2.+7.*r1*r2+r2**2.))/(r12-r1-r2)**7. &
      +(r12**2.+7.*r12*(r1+r2)+6*(r1**2.+7.*r1*r2+r2**2.))/(r12+r1+r2)**7. &
      -(r12**2.+7.*r12*(r1-r2)+6*(r1**2.-7.*r1*r2+r2**2.))/(r12+r1-r2)**7. &
      -(r12**2.-7.*r12*(r1-r2)+6*(r1**2.-7.*r1*r2+r2**2.))/(r12-r1+r2)**7. )

return
!------------------------------------------------------------------------------------------------------!
end subroutine hamaker_sphere_sphere



subroutine hamaker_plate_plate(h_12, A1, A2, Uatt, Urep)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: h_12, A1, A2
real(8), intent(out) :: Urep, Uatt
real(8)              :: hamaker_constant = 0.d0
real(8)              :: temp_place_holder = 0.d0
!------------------------------------------------------------------------------------------------------!
Urep = 0.d0
Uatt = 0.d0

temp_place_holder = h_12

hamaker_constant = SQRT(A1*A2)

Uatt = -hamaker_constant
return
!------------------------------------------------------------------------------------------------------!
end subroutine hamaker_plate_plate

!--------------------------------------------------------------------!
end module force_fields
