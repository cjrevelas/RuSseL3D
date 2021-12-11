!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_delta(numnp, qmx_interp_mg, ns_gr_conv, num_gpoints, gpid, delta_numer, gp_init_value, volnp)
!------------------------------------------------------------------------------------------------------!
use geometry,     only: box_lo, box_hi, xc, ndm
use iofiles,      only: gp_filename
use write_helper, only: adjl
use constants,    only: m3_to_A3
use error_handing
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: numnp, num_gpoints, ns_gr_conv
integer, intent(in), dimension(num_gpoints) :: gpid
integer                                     :: ii, iog, jj

real(8), intent(in), dimension(num_gpoints)        :: gp_init_value, delta_numer
real(8), intent(in), dimension(numnp)              :: volnp
real(8), intent(in), dimension(ns_gr_conv+1,numnp) :: qmx_interp_mg
!------------------------------------------------------------------------------------------------------!
! Update gnodes.in.lammpstrj input file
iog = 19
open(unit=iog, file=gp_filename)
write(iog,'(A14)') "ITEM: TIMESTEP"
write(iog,'(A1)')  '0'
write(iog,'(A21)') "ITEM: NUMBER OF ATOMS"
write(iog,'(I10)')  num_gpoints
write(iog,'(A25)') "ITEM: BOX BOUNDS PP PP PP"
do ii = 1, 3
    write(iog,'(2(F20.9))') box_lo(ii), box_hi(ii)
enddo
write(iog,'(A28)') "ITEM: ATOMS id type xu yu zu"

do ii = 1, num_gpoints
    write(iog,'(I10,2(1X,E20.9),3(1X,F20.9))') gpid(ii), gp_init_value(ii), delta_numer(ii), (xc(jj,gpid(ii)),jj=1,ndm)
enddo
close(iog)

open(unit=5, file = delta_out)
write(5,'(4(2X,A16))') "gpid", "qmx(rgi,Ng)", "Numerical delta", "Analytic delta"
do ii = 1, num_gpoints
    write(5,'(2X,I16,3(2X,E16.9))') gpid(ii), qmx_interp_mg(ns_gr_conv+1,gpid(ii)), delta_numer(ii), 1.d0 / volnp(gpid(ii)) * m3_to_A3
enddo
close(5)
!------------------------------------------------------------------------------------------------------!
end subroutine export_delta
