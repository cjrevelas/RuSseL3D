subroutine convolution(numnp, chainlen, ns, koeff, q1, q2, phia)
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in) :: numnp, ns
integer             :: k1, time_step

real(8), intent(in)                        :: chainlen
real(8), intent(in), dimension(ns+1)       :: koeff
real(8), intent(in), dimension(ns+1,numnp) :: q1, q2
real(8), intent(out), dimension(numnp)     :: phia
real(8)                                    :: summer
!-------------------------------------------------------------------------------------------!
do k1 = 1, numnp
    summer = 0.d0

    do time_step = 1, ns+1
        summer = summer + koeff(time_step) * q1(time_step,k1) * q2(ns+2-time_step,k1)
    enddo

    phia(k1) = summer / chainlen
enddo
!-------------------------------------------------------------------------------------------!
end subroutine convolution
