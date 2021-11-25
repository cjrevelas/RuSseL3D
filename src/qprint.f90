subroutine qprint
!--------------------------------------------------------------------!
use xdata
use mdata
use constants
!--------------------------------------------------------------------!   
implicit none
!--------------------------------------------------------------------!
integer :: time_step
real(8) :: iqf_final
!--------------------------------------------------------------------!	

! APS. Old section
!character(len=20) :: frmt
!write(frmt,'("(3E20.9,",I4,"E20.9)")') ns+1
!open(unit=363, file = 'qfree.out.txt')
!do i1 = 1, numnp
!     write (363,frmt) (xc(ii2,i1), ii2 = 1, ndm), (qf_final(i1,time_step), time_step = 1, ns+1)
!enddo
!close(363)
! / Old section

open(unit=363, file = 'qfree.out.txt')
do i1 = 1, numnp
    do ii2 = 1, ndm
        write(363,'(E20.9)',advance='no') xc(ii2,i1)
    enddo
    do time_step = 1, ns+1
        iqf_final = qf_final(i1,time_step)
        if (dabs(iqf_final)<tol) then
            iqf_final = 0.d0
        endif
        write(363,'(E20.9)',advance='no') iqf_final
    enddo
    write(363,*)
enddo
close(363)

return
!--------------------------------------------------------------------!
end subroutine qprint 
