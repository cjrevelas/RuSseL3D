!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module error_handling_mod
!----------------------------------------------------------------------------------------------------------------------------!
use iofiles_mod
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
logical :: FILE_EXISTS
integer :: ERROR_TYPE
character(100) :: ERROR_MESSAGE
!----------------------------------------------------------------------------------------------------------------------------!
  contains
    subroutine ExitWithError(STOP_SIGNAL, ERROR_TYPE, SCREEN, ERROR_MESSAGE)
      implicit none
      character(100) :: ERROR_MESSAGE
      integer        :: ERROR_TYPE
      integer        :: ioe = 354
      integer        :: STOP_SIGNAL
      integer        :: SCREEN

      open(unit=ioe, file = IO_errorFile, position = "append")

      if (ERROR_TYPE.eq.1) then
        write(ioe,'(A12)',advance='no') "PROBLEM BETWEEN COMPUTER AND SCREEN:"
      elseif (ERROR_TYPE.eq.2) then
        write(ioe,'(A12)',advance='no') "RUN ERROR:"
      endif

      write(ioe,*) ERROR_MESSAGE

      if (SCREEN.eq.1) then
        if (ERROR_TYPE.eq.1) then
          write(*,'(A33)',advance='no') "PROBLEM BETWEEN SCREEN AND CHAIR:"
        elseif (ERROR_TYPE.eq.2) then
          write(*,'(A12)',advance='no') "RUN ERROR:"
        endif
        write(*,*) ERROR_MESSAGE
      endif

      close(ioe)

      if (STOP_SIGNAL.eq.1) then
        write(*,*) " Exiting.."
        STOP
      endif
    end subroutine ExitWithError
!----------------------------------------------------------------------------------------------------------------------------!
end module error_handling_mod
