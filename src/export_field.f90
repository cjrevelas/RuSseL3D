subroutine export_field(wa, wa_new, wa_mix)
!-----------------------------------------------------------------------------------------------------------!
use geometry,     only: numnp, xc
use iofiles,      only: field_profile
use write_helper, only: adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: kk

real(8), intent(in), dimension(numnp) :: wa, wa_new, wa_mix
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting field.",40)
open (unit=120, file = field_profile)
write(120,'(7(A19))') "np", "x", "y", "z", "wa", "wa_new", "wa_mix"
do kk = 1, numnp
   write(120,'(I19,6(E19.9E3))') kk, xc(1,kk), xc(2,kk), xc(3,kk), wa(kk), wa_new(kk), wa_mix(kk)
enddo
close(120)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_field
