integer function get_sys_time()
!----------------------------------------------------------------------------------------------------------!
use iofiles, only: time
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer(4) :: status
!----------------------------------------------------------------------------------------------------------!
status = system("date +%s > d.time.out.txt")
open(unit=78, file = time)
read(78, '(I20)')get_sys_time
close(78)
status = system("rm d.time.out.txt")
return
!----------------------------------------------------------------------------------------------------------!
end function get_sys_time
