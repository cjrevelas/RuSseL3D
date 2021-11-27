subroutine convolution(q_final, phia)
!-------------------------------------------------------------------------------------------! 
use xdata
!-------------------------------------------------------------------------------------------! 
implicit none
!-------------------------------------------------------------------------------------------! 
integer :: k1, time_step                               !CJR

real(8), intent(in), dimension(numnp,ns+1) :: q_final  !CJR
real(8), intent(out), dimension(numnp)     :: phia     !CJR
!-------------------------------------------------------------------------------------------! 
do k1 = 1, numnp
    sum = 0.d0

    do time_step = 1, ns+1
        sum = sum + koeff(time_step)              *q_final(k1,time_step)*q_final(k1,ns+2-time_step)   !CJR
       !sum = sum + koeff(time_step)*ds(time_step)*q_final(k1,time_step)*q_final(k1,ns+2-time_step)   !CJR
    enddo

    phia(k1) = sum   !CJR
enddo


!-------------------------------------------------------------------------------------------! 
end subroutine convolution 
