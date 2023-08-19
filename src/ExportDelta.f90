!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportDelta(qmx_interp_mg, deltaNumerical, graftPointValue)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: numNodes, boxLow, boxHigh, nodeCoord, numDimensions
use parser_vars_mod,  only: numConvolPointsGrafted
use delta_mod,        only: targetNumGraftedChains, graftPointId
use arrays_mod,       only: nodeVolume
use iofiles_mod,      only: IO_graftFile
use write_helper_mod, only: adjl
use constants_mod,    only: m3_to_A3
use error_handling_mod
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer             :: ii, iog, jj

real(8), intent(in), dimension(targetNumGraftedChains)            :: graftPointValue, deltaNumerical
real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qmx_interp_mg
!------------------------------------------------------------------------------------------------------!
! Update gnodes.in.lammpstrj input file
iog = 19
open(unit=iog, file=IO_graftFile)
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
  write(iog,'(I10,2(1X,E20.9),3(1X,F20.9))') graftPointId(ii), graftPointValue(ii), deltaNumerical(ii), &
                                             (nodeCoord(jj,graftPointId(ii)), jj=1, numDimensions)
enddo
close(iog)

open(unit=5, file = IO_delta)  ! CJR: how does this work?? IO_delta has not been invoked from iofiles_mod!
write(5,'(4(2X,A16))') "gpid", "qmx(rgi,Ng)", "Numerical delta", "Analytic delta"
do ii = 1, targetNumGraftedChains
  write(5,'(2X,I16,3(2X,E16.9))') graftPointId(ii), qmx_interp_mg(numConvolPointsGrafted+1,graftPointId(ii)), &
                                  deltaNumerical(ii), 1.0d0 / nodeVolume(graftPointId(ii)) * m3_to_A3
enddo
close(5)
!------------------------------------------------------------------------------------------------------!
end subroutine ExportDelta
