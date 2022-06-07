!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_contour_coeffs(ds, ns, coeff)
!------------------------------------------------------------------------------!
use iofiles_mod, only: IO_contourCoeffs
!------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: nn

real(8), intent(in), dimension(ns+1)  :: ds
real(8), intent(out), dimension(ns+1) :: coeff
real(8), dimension(ns+1)              :: xx
!------------------------------------------------------------------------------!
xx(1) = 0.0d0
do nn = 2, ns+1
  xx(nn) = xx(nn-1) + ds(nn)
enddo

coeff = 0.0d0
do nn = 2, ns, 2
  coeff(nn-1) = coeff(nn-1) + (xx(nn+1) - xx(nn-1))*(2.0d0*xx(nn-1) + xx(nn+1) - 3.0d0*xx(nn)) / (6.0d0*(xx(nn-1) - xx(nn)))

  coeff(nn)   = coeff(nn)   + (xx(nn-1) - xx(nn+1))**3.0d0 / (6.0d0*(xx(nn) - xx(nn-1))*(xx(nn) - xx(nn+1)))

  coeff(nn+1) = coeff(nn+1) + (xx(nn+1) - xx(nn-1))*(2.0d0*xx(nn+1) + xx(nn-1) - 3.0d0*xx(nn)) / (6.0d0*(xx(nn+1) - xx(nn)))
enddo

#ifdef DEBUG_OUTPUTS
open(unit=400, file = IO_contourCoeffs, position = 'append')
write(400,'(5(A17))')  "n", "s", "ds", "coeff", "coeff_reduced"
do nn = 1, ns+1
  write(400,'(I17, 4(E17.9))') nn, xx(nn), ds(nn), coeff(nn), coeff(nn)/ds(nn)
enddo

write(400,*)
write(400,*) "-----------------------------------------------------------------------------------"
write(400,*) "-----------------------------------------------------------------------------------"
write(400,*)
close(400)
#endif

return
!------------------------------------------------------------------------------!
end subroutine compute_contour_coeffs
