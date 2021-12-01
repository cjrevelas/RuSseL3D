subroutine part_fun(numnp, ns, q_final, part_func)
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in)                        :: numnp, ns
integer                                    :: k1

real(8), intent(in), dimension(numnp,ns+1) :: q_final
real(8), intent(out)                       :: part_func
real(8)                                    :: sum_f, Q
real(8), dimension(numnp)                  :: q_last
!-------------------------------------------------------------------------------------------! 
sum_f = 0.d0

do k1 = 1, numnp
    q_last(k1) = q_final(k1,ns+1)
enddo

call spat_3d(q_last, sum_f, Q)

part_func = Q

return
!-------------------------------------------------------------------------------------------!
end subroutine part_fun
