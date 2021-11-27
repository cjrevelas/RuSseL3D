subroutine convolution(q1, q2, phia)
!-------------------------------------------------------------------------------------------! 
use xdata
!-------------------------------------------------------------------------------------------! 
implicit none
!-------------------------------------------------------------------------------------------! 
integer :: k1, time_step

real(8), intent(in), dimension(numnp,ns+1) :: q1, q2
real(8), intent(out), dimension(numnp)     :: phia
real(8)                                    :: summer
!-------------------------------------------------------------------------------------------! 
do k1 = 1, numnp
    summer = 0.d0

    do time_step = 1, ns+1
        summer = summer + koeff(time_step)              *q1(k1,time_step)*q2(k1,ns+2-time_step)
       !summer = summer + koeff(time_step)*ds(time_step)*q1(k1,time_step)*q2(k1,ns+2-time_step)
    enddo

    phia(k1) = summer
enddo


!-------------------------------------------------------------------------------------------! 
end subroutine convolution 
