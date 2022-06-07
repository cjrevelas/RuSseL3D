!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module write_helper_mod
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
integer, parameter :: PLOG    = 0
integer, parameter :: PSCREEN = 1
integer, parameter :: PBOTH   = 2

character(len=100) :: MESSAGE
!----------------------------------------------------------------------------------------------------------------------------!
  contains
    subroutine output(TARGET, iow, MESSAGE)
      implicit none
      character(len=100) :: MESSAGE
      integer            :: TARGET
      integer            :: iow

      if (TARGET.eq.0) then
        write(iow,*) MESSAGE
      elseif (TARGET.eq.1) then
        write(6,*) MESSAGE
      elseif (TARGET.eq.2) then
        write(iow,*) MESSAGE
        write(6,*) MESSAGE
      else
        write(6,*) "WRITE_HELPER: WRONG TARGET TO output FUNCTION"
        write(6,*) "    CHOOSE BETWEEN 0, 1 and 2"
      endif
    end subroutine output

    function adjl(string,length) result(r)
      implicit none
      character(len=*)      :: string
      integer               :: length
      character(len=length) :: r

      r = ADJUSTL(string)
    end function adjl

    logical function export(export_freq, iter, convergence)
      implicit none
      integer, intent(in) :: export_freq, iter
      logical, intent(in) :: convergence
      export=.false.
      if (export_freq>0) then
        if ((MOD(iter,export_freq).eq.0).or.convergence) export=.true.
      endif
    end function export
!----------------------------------------------------------------------------------------------------------------------------!
end module write_helper_mod
