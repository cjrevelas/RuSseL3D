integer function get_sys_time()
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer(4) :: status
!----------------------------------------------------------------------------------------------------------!
status = system("date +%s > time.temp")
open(unit=78, file = 'time.temp')
read(78, '(I20)')get_sys_time
close(78)
return
!----------------------------------------------------------------------------------------------------------!
end function get_sys_time