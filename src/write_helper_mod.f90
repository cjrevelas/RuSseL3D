module write_helper
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!

integer,parameter :: PLOG = 0
integer,parameter :: PSCREEN = 1
integer,parameter :: PBOTH = 2

CHARACTER(len=100) :: MESSAGE

Contains

subroutine output(TARGET, iow, MESSAGE)
   character(len=100) MESSAGE
   integer :: TARGET
   integer :: iow

   if (TARGET.eq.0) then
      write(iow,*)MESSAGE
   elseif (TARGET.eq.1) then
      write(6,*)MESSAGE
   elseif (TARGET.eq.2) then
      write(iow,*)MESSAGE
      write(6,*)MESSAGE
   else
      write(6,*)"WRITE_HELPER: WRONG TARGET TO output FUNCTION"
      write(6,*)"    CHOOSE BETWEEN 0, 1 and 2"
   endif


end subroutine output

function adjl(string,length) result(r)
   character(len=*) :: string
   integer          :: length
   character(len=length) :: r
   r = adjustl(string)
end function adjl

!--------------------------------------------------------------------!
end module write_helper
