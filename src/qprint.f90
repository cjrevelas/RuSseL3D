subroutine qprint
!--------------------------------------------------------------------!
use xdata
use mdata
!--------------------------------------------------------------------!   
implicit none
!--------------------------------------------------------------------!
character(len=20) :: frmt
!--------------------------------------------------------------------!	
write(frmt,'("(3E20.9,",I4,"E20.9)")') ns+1

open(unit=363, file = 'qfree.txt')

do i1 = 1, numnp
     write (363,frmt) (xc(ii2,i1), ii2 = 1, ndm), (qf_final(i1,time_step), time_step = 1, ns+1)
enddo  

close(363)

return 
!--------------------------------------------------------------------!
end subroutine qprint 
