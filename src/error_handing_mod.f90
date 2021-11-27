module error_handing
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
character(100) :: ERROR_MESSAGE
integer        :: ERROR_TYPE
logical        :: FILE_EXISTS
!--------------------------------------------------------------------!
contains

subroutine exit_with_error(STOP_SIGNAL, ERROR_TYPE, SCREEN, ERROR_MESSAGE)
    character(100) :: ERROR_MESSAGE
    integer        :: ERROR_TYPE
    integer        :: ioe = 354
    integer        :: STOP_SIGNAL
    integer        :: SCREEN


    open(unit=ioe, file = 'error.out.txt', position = 'append')

    if (ERROR_TYPE.eq.1) then
        write(ioe,'(A12)',advance='no')'PROBLEM BETWEEN COMPUTER AND SCREEN:'
    elseif (ERROR_TYPE.eq.2) then
        write(ioe,'(A12)',advance='no')'RUN ERROR:'
    endif

    write(ioe,*)ERROR_MESSAGE

    if (SCREEN.eq.1) then
        if (ERROR_TYPE.eq.1) then
            write(*,'(A33)',advance='no')'PROBLEM BETWEEN SCREEN AND CHAIR:'
        elseif (ERROR_TYPE.eq.2) then
            write(*,'(A12)',advance='no')'RUN ERROR:'
        endif
        write(*,*)ERROR_MESSAGE
    endif

    close(ioe)

    if (STOP_SIGNAL.eq.1) then
        write(*,*)' Exiting..'
        STOP
    endif
end subroutine
!--------------------------------------------------------------------!
end module error_handing
