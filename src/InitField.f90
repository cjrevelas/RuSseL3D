!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine InitField(Ufield, ww)
!------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: polymerAlpha, polymerSigma, fieldInitScheme, kapa, numNanopFaces, &
                            molarBulkDensity, beta, wallDistance, nanopCenter, plateSigma,    &
                            plateAlpha, nanopRadiusEff, nanopSigma, nanopAlpha
use geometry_mod,     only: numNodes, isDirichletFace, boxLow, boxHigh, nodeCoord, nodeBelongsToDirichletFace
use error_handling_mod
use write_helper_mod, only: adjl
use force_fields_mod, only: hamaker_sphere_plate, hamaker_sphere_sphere
use iofiles_mod,      only: IO_solidPotential, IO_fieldFile
use constants_mod,    only: n_avog, pi, m_to_A
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: kk, mm, nn

real(8), intent(out), dimension(numNodes) :: Ufield, ww
real(8)                                   :: number_density=0.0d0
real(8)                                   :: r_centers=0.0d0, r_center_surf=0.0d0, r_surf=0.0d0
real(8)                                   :: nanopRadiusActual=0.0d0, radius_pol=0.0d0
real(8)                                   :: Urep=0.0d0, Uatt=0.0d0
!------------------------------------------------------------------------------------------------------!
open(unit=211, file = IO_solidPotential)
write(211,'(5(2X,A16))') "r_center_surf", "r_surf", "Uatt", "Urep", "Utot"

! TODO: compute and add solid-solid (Hamaker) interactions
number_density = molarBulkDensity * n_avog
radius_pol     = (3.0d0 / (4.0d0 * pi * number_density))**(1.0d0/3.0d0) * m_to_A

Ufield = 0.0d0

do kk = 1, numNodes
  ! Loop over all dirichlet faces
  do mm = 1, 3
    do nn = 1, 2
      if (isDirichletFace(mm,nn)) then
        number_density = molarBulkDensity * n_avog
        radius_pol     = (3.0d0 / 4.0d0 / pi / number_density)**(1.0d0/3.0d0) * m_to_A

        if (nn.eq.1) then
          r_center_surf = nodeCoord(mm,kk) - boxLow(mm) + wallDistance
        elseif (nn.eq.2) then
          r_center_surf = boxHigh(mm) - nodeCoord(mm,kk) + wallDistance
        endif

        r_surf = r_center_surf - radius_pol

        call hamaker_sphere_plate(r_surf, radius_pol, polymerSigma, plateSigma(1), polymerAlpha, plateAlpha(1), Urep, Uatt)

        Urep = Urep*beta
        Uatt = Uatt*beta

        Ufield(kk) = Ufield(kk) + Urep + Uatt

        if ((Ufield(kk)-Ufield(kk)).gt.1.0d-8) then
          write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". NaN was changed to ",E16.9)') r_surf, Ufield(kk)
          call exitWithError(1,2,1,ERROR_MESSAGE)
        endif

        if (r_surf.lt.0.0d0) then
          write(ERROR_MESSAGE,'("Hamaker distance smaller than zero! (",E16.9,").")') r_surf
          call exitWithError(1, 2, 1, ERROR_MESSAGE)
        endif

        write(211,'(5(2X,E16.9E2))') r_center_surf, r_surf, Uatt, Urep, Urep+Uatt
      endif
    enddo
  enddo

  ! Loop over all nanoparticle faces
  do mm = 1, numNanopFaces
    r_centers         = DSQRT((nodeCoord(1,kk)-nanopCenter(1,mm))**2.0d0 + (nodeCoord(2,kk)-nanopCenter(2,mm))**2.0d0 + (nodeCoord(3,kk)-nanopCenter(3,mm))**2.0d0)
    nanopRadiusActual = nanopRadiusEff(mm) - wallDistance
    r_surf            = r_centers - radius_pol - nanopRadiusActual

    call hamaker_sphere_sphere(r_surf, radius_pol, nanopRadiusActual, polymerSigma, nanopSigma(mm), polymerAlpha, nanopAlpha(mm), Urep, Uatt)

    Urep = Urep * beta
    Uatt = Uatt * beta

    Ufield(kk) = Ufield(kk) + Urep + Uatt

    if ((Ufield(kk)-Ufield(kk)).gt.1.0d-8) then
      write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". NaN was changed to ",E16.9)') r_surf, Ufield(kk)
      call exitWithError(1, 2, 1, ERROR_MESSAGE)
    endif

    write(211,'(5(2X,E16.9E2))') r_centers, r_surf, Uatt, Urep, Urep+Uatt
  enddo
enddo

close(211)

if (fieldInitScheme.eq.0) then
  ww = 0.0d0
elseif (fieldInitScheme.eq.1) then
  open(unit=655, file = IO_fieldFile, Form='unformatted')
  read(655) ww
  close(655)
elseif (fieldInitScheme.eq.2) then
  do kk = 1, numNodes
    if (nodeBelongsToDirichletFace(kk)) then
      ww(kk) = -kapa
    else
      ww(kk) = 0.0d0
    endif
  enddo
endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine InitField
