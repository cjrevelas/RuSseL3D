!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_part_func_mx(numnp, ns, q_final, part_func)
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in)                        :: numnp, ns
integer                                    :: kk

real(8), intent(in), dimension(ns+1,numnp) :: q_final
real(8), intent(out)                       :: part_func
real(8)                                    :: sum_f, vol
real(8), dimension(numnp)                  :: q_last
!-------------------------------------------------------------------------------------------!
sum_f = 0.0d0

do kk = 1, numnp
    q_last(kk) = q_final(ns+1,kk)
enddo

call fem_integration(q_last, sum_f, part_func, vol)

return
!-------------------------------------------------------------------------------------------!
end subroutine compute_part_func_mx
