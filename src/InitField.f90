!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine InitField(uuField, ww)
!------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: polymerAlpha, polymerSigma, fieldInitScheme, kapa, numNanopFaces, &
                            molarBulkDensity, beta, wallDistance, nanopCenter, plateSigma,    &
                            plateAlpha, nanopRadiusEff, nanopSigma, nanopAlpha
use geometry_mod,     only: numNodes, isDirichletFace, boxLow, boxHigh, nodeCoord, nodeBelongsToDirichletFace
use error_handling_mod
use write_helper_mod, only: adjl
use force_fields_mod, only: HamakerSpherePlate, HamakerSphereSphere
use iofiles_mod,      only: IO_solidPotential, IO_fieldFile
use constants_mod,    only: n_avog, pi, m_to_A
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: node, axis, face

real(8), intent(out), dimension(numNodes) :: uuField, ww

real(8) :: numDensity=0.0d0
real(8) :: centerToCenterDistance=0.0d0, centerToSurfaceDistance=0.0d0, surfaceToSurfaceDistance=0.0d0
real(8) :: nanopRadiusActual=0.0d0, segmentRadius=0.0d0
real(8) :: uuRepulsive=0.0d0, uuAttractive=0.0d0
!------------------------------------------------------------------------------------------------------!
open(unit=211, file = IO_solidPotential)
write(211,'(5(2X,A16))') "centerToSurfaceDistance", "surfaceToSurfaceDistance", "attractive", "repulsive", "total"

numDensity    = molarBulkDensity * n_avog
segmentRadius = (3.0d0 / (4.0d0 * pi * numDensity))**(1.0d0/3.0d0) * m_to_A

uuField = 0.0d0

do node = 1, numNodes
  ! Loop over all Dirichlet faces
  do axis = 1, 3
    do face = 1, 2
      if (isDirichletFace(axis,face)) then
        numDensity    = molarBulkDensity * n_avog
        segmentRadius = (3.0d0 / 4.0d0 / pi / numDensity)**(1.0d0/3.0d0) * m_to_A

        if (face.eq.1) then
          centerToSurfaceDistance = nodeCoord(axis,node) - boxLow(axis) + wallDistance
        elseif (face.eq.2) then
          centerToSurfaceDistance = boxHigh(axis) - nodeCoord(axis,node) + wallDistance
        endif

        surfaceToSurfaceDistance = centerToSurfaceDistance - segmentRadius

        CALL HamakerSpherePlate(surfaceToSurfaceDistance, segmentRadius, polymerSigma, plateSigma(1), polymerAlpha, plateAlpha(1), uuRepulsive, uuAttractive)

        uuRepulsive  = uuRepulsive * beta
        uuAttractive = uuAttractive * beta

        uuField(node) = uuField(node) + uuRepulsive + uuAttractive

        if ((uuField(node)-uuField(node)).gt.1.0d-8) then
          write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". NaN was changed to ",E16.9)') surfaceToSurfaceDistance, uuField(node)
          CALL ExitWithError(1,2,1,ERROR_MESSAGE)
        endif

        if (surfaceToSurfaceDistance.lt.0.0d0) then
          write(ERROR_MESSAGE,'("Hamaker distance smaller than zero! (",E16.9,").")') surfaceToSurfaceDistance
          CALL ExitWithError(1, 2, 1, ERROR_MESSAGE)
        endif

        write(211,'(5(2X,E16.9E2))') centerToSurfaceDistance, surfaceToSurfaceDistance, uuAttractive, uuRepulsive, uuRepulsive + uuAttractive
      endif
    enddo
  enddo

  ! Loop over all nanoparticle faces
  do face = 1, numNanopFaces
    centerToCenterDistance   = DSQRT((nodeCoord(1,node)-nanopCenter(1,face))**2.0d0 + &
                                     (nodeCoord(2,node)-nanopCenter(2,face))**2.0d0 + &
                                     (nodeCoord(3,node)-nanopCenter(3,face))**2.0d0)

    nanopRadiusActual        = nanopRadiusEff(face)   - wallDistance
    surfaceToSurfaceDistance = centerToCenterDistance - segmentRadius - nanopRadiusActual

    CALL HamakerSphereSphere(surfaceToSurfaceDistance, segmentRadius, nanopRadiusActual, polymerSigma, nanopSigma(face), polymerAlpha, nanopAlpha(face), uuRepulsive, uuAttractive)

    uuRepulsive  = uuRepulsive  * beta
    uuAttractive = uuAttractive * beta

    uuField(node) = uuField(node) + uuRepulsive + uuAttractive

    if ((uuField(node)-uuField(node)).gt.1.0d-8) then
      write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". NaN was changed to ",E16.9)') surfaceToSurfaceDistance, uuField(node)
      CALL ExitWithError(1, 2, 1, ERROR_MESSAGE)
    endif

    write(211,'(5(2X,E16.9E2))') centerToCenterDistance, surfaceToSurfaceDistance, uuAttractive, uuRepulsive, uuRepulsive + uuAttractive
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
  do node = 1, numNodes
    if (nodeBelongsToDirichletFace(node)) then
      ww(node) = -kapa
    else
      ww(node) = 0.0d0
    endif
  enddo
endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine InitField
