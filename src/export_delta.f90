subroutine export_delta(numnp, qmx_interp_mg, ns_gr_conv, num_gpoints, gpid, delta_numer, gp_init_value, volnp)
!------------------------------------------------------------------------------------------------------!
use geometry,     only: box_lo, box_hi
use iofiles,      only: gp_filename
use write_helper, only: adjl
use error_handing
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: numnp, num_gpoints, ns_gr_conv
integer, intent(in), dimension(num_gpoints) :: gpid
integer                                     :: ii, iog

real(8), intent(in), dimension(num_gpoints)        :: gp_init_value
real(8), intent(in), dimension(numnp)              :: volnp
real(8), intent(in), dimension(ns_gr_conv+1,numnp) :: qmx_interp_mg
real(8), intent(out), dimension(num_gpoints)       :: delta_numer
!------------------------------------------------------------------------------------------------------!
!update gnodes.in.lammpstrj input file
iog = 19
open(unit=iog, file=gp_filename)
write(iog,'(A14)') "ITEM: TIMESTEP"
write(iog,'(A1)')  '0'
write(iog,'(A21)') "ITEM: NUMBER OF ATOMS"
write(iog,'(I10)')  num_gpoints
write(iog,'(A25)') "ITEM: BOX BOUNDS PP PP PP"
do ii = 1, 3
    write(iog,'(2(F19.15))') box_lo(ii), box_hi(ii)
enddo
write(iog,'(A28)') "ITEM: ATOMS id type xu yu zu"

do ii = 1, num_gpoints
    write(iog,'(I5,2(1X,E16.9))') gpid(ii), gp_init_value(ii), delta_numer(ii)
enddo
close(iog)

open(unit=5, file = delta_out)
write(5,'(A5,A16,A20,A17)') "gpid", "qmx(rgi,Ng)", "Numerical delta", "Analytic delta"
do ii = 1, num_gpoints
    write(5,'(I5,3(1X,E18.9))') gpid(ii), qmx_interp_mg(ns_gr_conv+1,gpid(ii)), delta_numer(ii), 1.d0 / volnp(gpid(ii)) * 1.e+30
enddo
close(5)
!------------------------------------------------------------------------------------------------------!
end subroutine export_delta
