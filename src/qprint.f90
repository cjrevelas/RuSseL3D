subroutine qprint(ns, q_final, q_type)
!--------------------------------------------------------------------!
use geometry
use constants
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: time_step, i1, ii2

character(4)                               :: q_type
character(20)                              :: file_name

real(8), intent(in), dimension(numnp,ns+1) :: q_final
real(8)                                    :: iq_final
!--------------------------------------------------------------------!
write(file_name,'(''q'',A4,''.out.txt'')') q_type

open(unit=363, file = file_name)
do i1 = 1, numnp
    do ii2 = 1, ndm
        write(363,'(E20.9)',advance='no') xc(ii2,i1)
    enddo
    do time_step = 1, ns+1
        iq_final = q_final(i1,time_step)
        if (dabs(iq_final)<tol) then
            iq_final = 0.d0
        endif
        write(363,'(E20.9)',advance='no') iq_final
    enddo
    write(363,*)
enddo
close(363)

return
!--------------------------------------------------------------------!
end subroutine qprint
