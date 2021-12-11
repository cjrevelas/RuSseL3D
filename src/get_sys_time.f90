!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

integer function get_sys_time()
!----------------------------------------------------------------------------------------------------------!
use iofiles_mod, only: time
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer(4) :: status
!----------------------------------------------------------------------------------------------------------!
status = system("date +%s > o.time")
open(unit=78, file = time)
read(78, '(I20)')get_sys_time
close(78)
status = system("rm o.time")
return
!----------------------------------------------------------------------------------------------------------!
end function get_sys_time
