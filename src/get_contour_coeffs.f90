!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine get_contour_coeffs(ds, ns, coeff)
!------------------------------------------------------------------------------!
use iofiles, only: contour_coeffs
!------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: nn

real(8), intent(in), dimension(ns+1)  :: ds
real(8), intent(out), dimension(ns+1) :: coeff
real(8), dimension(ns+1)              :: x
!------------------------------------------------------------------------------!
x(1) = 0.d0
do nn = 2, ns+1
   x(nn) = x(nn-1) + ds(nn)
enddo

coeff = 0.d0
do nn = 2, ns, 2
   coeff(nn-1) = coeff(nn-1) + ( x(nn+1)-x(nn-1) )*( 2*x(nn-1)+x(nn+1)-3*x(nn) ) &
      &                    / ( 6.d0 * ( x(nn-1) - x(nn) ))

   coeff(nn)   = coeff(nn)   + ( x(nn-1) - x(nn+1) )**3.d0 &
      &                    / ( 6.d0 * ( x(nn) - x(nn-1) ) * (x(nn) - x(nn+1) ) )

   coeff(nn+1) = coeff(nn+1) + ( x(nn+1)-x(nn-1) )*( 2*x(nn+1)+x(nn-1)-3*x(nn) ) &
      &                    / ( 6.d0 * ( x(nn+1) - x(nn) ))
enddo

#ifdef DEBUG_OUTPUTS
open(unit=400, file = contour_coeffs, position = 'append')
write(400,'(5(A17))')  "n", "s", "ds", "coeff", "coeff_reduced"
do nn = 1, ns+1
    write(400,'(I17, 4(E17.9))') nn, x(nn), ds(nn), coeff(nn), coeff(nn)/ds(nn)
enddo

write(400,*)
write(400,*) "-----------------------------------------------------------------------------------"
write(400,*) "-----------------------------------------------------------------------------------"
write(400,*)
close(400)
#endif

return
!------------------------------------------------------------------------------!
end subroutine get_contour_coeffs
