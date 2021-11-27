subroutine quadinterp_koef(koeff, x, ds, ns)
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: n, ns
real(8), intent(in), dimension(ns+1) :: x, ds
real(8), intent(out), dimension(ns+1) :: koeff

koeff = 0.d0

do n = 2, ns, 2
   koeff(n-1) = koeff(n-1) + ( x(n+1)-x(n-1) )*( 2*x(n-1)+x(n+1)-3*x(n) ) &
      &                    / ( 6.d0 * ( x(n-1) - x(n) ))

   koeff(n)   = koeff(n)   + ( x(n-1) - x(n+1) )**3.d0 &
      &                    / ( 6.d0 * ( x(n) - x(n-1) ) * (x(n) - x(n+1) ) )

   koeff(n+1) = koeff(n+1) + ( x(n+1)-x(n-1) )*( 2*x(n+1)+x(n-1)-3*x(n) ) &
      &                    / ( 6.d0 * ( x(n+1) - x(n) ))
enddo

!#ifdef DEBUG_OUTPUTS
open(unit=400, file = 'quad_interp_coeffs.out.txt')
write(400,'(5(A17))')      "n", "s",  "ds", "coeff",  "coeff_reduced"
do n = 1, ns+1
    write(400,'(I17, 4(E17.9))')n, x(n), ds(n), koeff(n), koeff(n)/ds(n)
enddo
close(400)
!#endif

return
!--------------------------------------------------------------------!
end subroutine quadinterp_koef
