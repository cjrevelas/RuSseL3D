subroutine quadinterp_koef(ds, ns, koeff)
!------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: n

real(8), intent(in), dimension(ns+1)  :: ds
real(8), intent(out), dimension(ns+1) :: koeff
real(8), dimension(ns+1)  :: x

x(1) = 0.d0
do n = 2, ns+1
   x(n) = x(n-1) + ds(n)
enddo
!------------------------------------------------------------------------------!
koeff = 0.d0

do n = 2, ns, 2
   koeff(n-1) = koeff(n-1) + ( x(n+1)-x(n-1) )*( 2*x(n-1)+x(n+1)-3*x(n) ) &
      &                    / ( 6.d0 * ( x(n-1) - x(n) ))

   koeff(n)   = koeff(n)   + ( x(n-1) - x(n+1) )**3.d0 &
      &                    / ( 6.d0 * ( x(n) - x(n-1) ) * (x(n) - x(n+1) ) )

   koeff(n+1) = koeff(n+1) + ( x(n+1)-x(n-1) )*( 2*x(n+1)+x(n-1)-3*x(n) ) &
      &                    / ( 6.d0 * ( x(n+1) - x(n) ))
enddo

#ifdef DEBUG_OUTPUTS
open(unit=400, file = 'quad_interp_coeffs.out.txt', position = 'append')
write(400,'(5(A17))')      "n", "s",  "ds", "coeff",  "coeff_reduced"
do n = 1, ns+1
    write(400,'(I17, 4(E17.9))')n, x(n), ds(n), koeff(n), koeff(n)/ds(n)
enddo

write(400,*)
write(400,*) "-----------------------------------------------------------------------------------"
write(400,*) "-----------------------------------------------------------------------------------"
write(400,*)
close(400)
#endif

return
!------------------------------------------------------------------------------!
end subroutine quadinterp_koef
