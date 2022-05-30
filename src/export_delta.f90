!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_delta(numNodes, qmx_interp_mg, numConvolPointsGrafted, targetNumGraftedChains, &
                        gpid, delta_numer, gp_init_value, volnp)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: boxLow, boxHigh, nodeCoord, numDimensions
use iofiles_mod,      only: graftFile
use write_helper_mod, only: adjl
use constants_mod,    only: m3_to_A3
use error_handing_mod
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                                    :: numNodes, targetNumGraftedChains, numConvolPointsGrafted
integer, intent(in), dimension(targetNumGraftedChains) :: gpid
integer                                                :: ii, iog, jj

real(8), intent(in), dimension(targetNumGraftedChains)            :: gp_init_value, delta_numer
real(8), intent(in), dimension(numNodes)                          :: volnp
real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qmx_interp_mg
!------------------------------------------------------------------------------------------------------!
! Update gnodes.in.lammpstrj input file
iog = 19
open(unit=iog, file=graftFile)
write(iog,'(A14)') "ITEM: TIMESTEP"
write(iog,'(A1)')  '0'
write(iog,'(A21)') "ITEM: NUMBER OF ATOMS"
write(iog,'(I10)')  targetNumGraftedChains
write(iog,'(A25)') "ITEM: BOX BOUNDS PP PP PP"
do ii = 1, 3
  write(iog,'(2(F20.9))') boxLow(ii), boxHigh(ii)
enddo
write(iog,'(A28)') "ITEM: ATOMS id type xu yu zu"

do ii = 1, targetNumGraftedChains
  write(iog,'(I10,2(1X,E20.9),3(1X,F20.9))') gpid(ii), gp_init_value(ii), delta_numer(ii), (nodeCoord(jj,gpid(ii)),jj=1,numDimensions)
enddo
close(iog)

open(unit=5, file = delta_out)
write(5,'(4(2X,A16))') "gpid", "qmx(rgi,Ng)", "Numerical delta", "Analytic delta"
do ii = 1, targetNumGraftedChains
  write(5,'(2X,I16,3(2X,E16.9))') gpid(ii), qmx_interp_mg(numConvolPointsGrafted+1,gpid(ii)), delta_numer(ii), 1.0d0 / volnp(gpid(ii)) * m3_to_A3
enddo
close(5)
!------------------------------------------------------------------------------------------------------!
end subroutine export_delta
