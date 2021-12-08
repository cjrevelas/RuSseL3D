subroutine get_part_func(numnp, ns, q_final, part_func)
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in)                        :: numnp, ns
integer                                    :: k1

real(8), intent(in), dimension(ns+1,numnp) :: q_final
real(8), intent(out)                       :: part_func
real(8)                                    :: sum_f, Q, vol
real(8), dimension(numnp)                  :: q_last
!-------------------------------------------------------------------------------------------!
sum_f = 0.d0

do k1 = 1, numnp
    q_last(k1) = q_final(ns+1,k1)
enddo

call spat_3d(q_last, sum_f, Q, vol)

part_func = Q

return
!-------------------------------------------------------------------------------------------!
end subroutine get_part_func
