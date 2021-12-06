subroutine simpsonkoef_S(ds_ave, ns, koeff)
!------------------------------------------------------------------------------------!
use iofiles
!------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: n

real(8), intent(in)                   :: ds_ave
real(8), intent(out), dimension(ns+1) :: koeff
!------------------------------------------------------------------------------------!
!*************************************************************!
! For constant step size the algorithm reduces to Simpson 1/3 !
! rule depicted in the code snipet below:                     !
!                                                             !
! s values:      0                                        1   !
!                                                             !
! Nodal points:  0    1   2   3   4   ...   ns-2   ns-1   ns  !
!                                                             !
! Coefficients: 1/3  4/3 2/3 4/3 2/3  ..     2/3   4/3    1/3 !
!                                                             !
!*************************************************************!
koeff(1)    = 1.d00/3.d00
koeff(ns+1) = 1.d00/3.d00
koeff(ns)   = 4.d00/3.d00

do n = 1, ns/2
     koeff(2*n)   = 4.d00/3.d00
     koeff(2*n+1) = 2.d00/3.d00
enddo

koeff(ns+1) = 1.d00/3.d00

do n = 1, ns+1
    koeff(n) = koeff(n) * ds_ave
enddo

open(unit=400, file = simpson)
write(400,'(5(A17))')   "n", "s", "ds", "coeff", "coeff_reduced"
do n = 1, ns+1
    write(400,'(I17, 4(E17.9))')n, (n-1)*ds_ave, ds_ave, koeff(n), koeff(n)/ds_ave
enddo
close(400)

return
!------------------------------------------------------------------------------------!
end subroutine simpsonkoef_S
