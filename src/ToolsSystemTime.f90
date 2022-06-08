!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

integer function ToolsSystemTime()
!----------------------------------------------------------------------------------------------------------!
use iofiles_mod, only: IO_time
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer(4) :: status
!----------------------------------------------------------------------------------------------------------!
status = system("date +%s > o.time")
open(unit=78, file = IO_time)
read(78, '(I20)') ToolsSystemTime
close(78)
status = system("rm o.time")
return
!----------------------------------------------------------------------------------------------------------!
end function ToolsSystemTime
