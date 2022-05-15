!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_field(Ufield, ww)
!------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: a_pol, field_init_scheme, kapa, numNanoparticleFaces,  &
&                       molarBulkDensity, sigma_pol, beta, wall_distance, center_np,        &
&                       sigma_plate, radius_np_eff, sigma_np, A_np, A_plate
use geometry_mod,     only: numnp, isDirichletFace, boxLow, boxHigh, xc, nodeBelongsToDirichletFace
use error_handing_mod
use write_helper_mod, only: adjl
use force_fields_mod, only: hamaker_sphere_plate, hamaker_sphere_sphere
use iofiles_mod,      only: usolid, fieldFile
use constants_mod,    only: n_avog, pi, m_to_A
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: kk, mm, nn

real(8), intent(out), dimension(numnp) :: Ufield, ww
real(8)                                :: number_density=0.0d0
real(8)                                :: r_centers=0.0d0, r_center_surf=0.0d0, r_surf=0.0d0
real(8)                                :: radius_np_actual=0.0d0, radius_pol=0.0d0
real(8)                                :: Urep=0.0d0, Uatt=0.0d0
!------------------------------------------------------------------------------------------------------!
open(unit=211, file = usolid)
write(211,'(5(2X,A16))') "r_center_surf", "r_surf", "Uatt", "Urep", "Utot"

! TODO: compute and add solid-solid (Hamaker) interactions

Ufield = 0.0d0

do kk = 1, numnp
  ! Loop over all dirichlet faces
  do mm = 1, 3
    do nn = 1, 2
      if (isDirichletFace(mm,nn)) then
        number_density = molarBulkDensity * n_avog
        radius_pol     = (3.0d0 / 4.0d0 / pi / number_density)**(1.0d0/3.0d0) * m_to_A

        if (nn.eq.1) then
          r_center_surf = xc(mm,kk) - boxLow(mm) + wall_distance
        elseif (nn.eq.2) then
          r_center_surf = boxHigh(mm) - xc(mm,kk) + wall_distance
        endif

        r_surf = r_center_surf - radius_pol

        call hamaker_sphere_plate(r_surf, radius_pol, sigma_pol, sigma_plate(1), A_pol, A_plate(1), Urep, Uatt)

        Urep = Urep*beta
        Uatt = Uatt*beta

        Ufield(kk) = Ufield(kk) + Urep + Uatt

        if ((Ufield(kk)-Ufield(kk)).gt.1.0d-8) then
          write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". &
                              & NaN was changed to ",E16.9)') r_surf, Ufield(kk)
          call exit_with_error(1,2,1,ERROR_MESSAGE)
        endif

        if (r_surf.lt.0.0d0) then
          write(ERROR_MESSAGE,'("Hamaker distance smaller than zero! (",E16.9,").")') r_surf
          call exit_with_error(1,2,1,ERROR_MESSAGE)
        endif

        write(211,'(5(2X,E16.9E2))') r_center_surf, r_surf, Uatt, Urep, Urep+Uatt
      endif
    enddo
  enddo

  ! Loop over all nanoparticle faces
  ! TODO: the following two lines can be moved outside of the outer loop
  number_density = molarBulkDensity * n_avog
  radius_pol     = (3.0d0 / 4.0d0 / pi / number_density)**(1.0d0/3.0d0) * m_to_A

  do mm = 1, numNanoparticleFaces
    r_centers        = DSQRT((xc(1,kk)-center_np(1,mm))**2.0d0 + (xc(2,kk)-center_np(2,mm))**2.0d0 + (xc(3,kk)-center_np(3,mm))**2.0d0)
    radius_np_actual = radius_np_eff(mm) - wall_distance
    r_surf           = r_centers - radius_pol - radius_np_actual

    call hamaker_sphere_sphere(r_surf, radius_pol, radius_np_actual, sigma_pol, sigma_np(mm), A_pol, A_np(mm), Urep, Uatt)

    Urep = Urep * beta
    Uatt = Uatt * beta

    Ufield(kk) = Ufield(kk) + Urep + Uatt

    if ((Ufield(kk)-Ufield(kk)).gt.1.0d-8) then
      write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". NaN was changed to ",E16.9)') r_surf, Ufield(kk)
      call exit_with_error(1,2,1,ERROR_MESSAGE)
    endif

    write(211,'(5(2X,E16.9E2))') r_centers, r_surf, Uatt, Urep, Urep+Uatt
  enddo
enddo

close(211)

if (field_init_scheme.eq.0) then
  ww = 0.0d0
elseif (field_init_scheme.eq.1) then
  open(unit=655, file = fieldFile, Form='unformatted')
  read(655) ww
  close(655)
elseif (field_init_scheme.eq.2) then
  do kk = 1, numnp
    if (nodeBelongsToDirichletFace(kk)) then
      ww(kk) = -kapa
    else
      ww(kk) = 0.0d0
    endif
  enddo
endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine init_field
