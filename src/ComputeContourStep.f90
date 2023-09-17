!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ComputeContourStep(dsAve, xsCrit, chainLength, nsEdw)
!------------------------------------------------------------------------------------------------------!
use constants_mod,    only: pi
use parser_vars_mod,  only: iow
use write_helper_mod, only: adjl
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(out) :: nsEdw
integer              :: nsPartOne, nsPartTwo

real(8), intent(in)  :: dsAve, xsCrit, chainLength
real(8)              :: dsMax
!------------------------------------------------------------------------------------------------------!
nsPartOne = 2 * NINT(0.5d0 * xsCrit / dsAve)

if (chainLength.ge.xsCrit) then
  dsMax     = (xsCrit * (1.0d0 - DCOS(pi * (DBLE(nsPartOne+1)-1.0d0) / (DBLE(nsPartOne) * 2.0d0)))) &
            - (xsCrit * (1.0d0 - DCOS(pi * (DBLE(nsPartOne)  -1.0d0) / (DBLE(nsPartOne) * 2.0d0))))
  nsPartTwo = 2 * NINT(0.5d0*(chainLength - xsCrit)/dsMax)
  nsEdw     = nsPartOne + nsPartTwo + 1

  write(iow,'(3X,A40)')adjl("Critical point smaller than chainlength.",40)
  write(6  ,'(3X,A40)')adjl("Critical point smaller than chainlength.",40)
  write(iow,'(3X,A40,F16.1)')adjl("Maximum contour step size:",40), dsMax
  write(6  ,'(3X,A40,F16.1)')adjl("Maximum contour step size:",40), dsMax
  write(iow,'(3X,A40,I16)')adjl("Contour steps before critical point:",40), nsPartOne
  write(6  ,'(3X,A40,I16)')adjl("Contour steps before critical point:",40), nsPartOne
  write(iow,'(3X,A40,I16)')adjl("Contour steps after critical point:",40), nsPartTwo
  write(6  ,'(3X,A40,I16)')adjl("Contour steps after critical point:",40), nsPartTwo
  write(iow,'(3X,A40,I16)')adjl("Total number of contour steps:",40), nsEdw
  write(6  ,'(3X,A40,I16)')adjl("Total number of contour steps:",40), nsEdw
else
  nsEdw = CEILING( 1.0d0 + 2.0d0 * xsCrit / (dsAve * pi) * ACOS(1.0d0 - chainLength / xsCrit) )
  dsMax = (xsCrit * (1.0d0 - DCOS(pi * (DBLE(nsEdw+1)-1.0d0) / (DBLE(nsPartOne) * 2.0d0)))) &
        - (xsCrit * (1.0d0 - DCOS(pi * (DBLE(nsEdw)  -1.0d0) / (DBLE(nsPartOne) * 2.0d0))))
  write(iow,'(3X,A40)')adjl("Critical point larger than chainlength.",40)
  write(6  ,'(3X,A40)')adjl("Critical point larger than chainlength.",40)
  write(iow,'(3X,A40,F16.1)')adjl("Maximum contour step size:",40), dsMax
  write(6  ,'(3X,A40,F16.1)')adjl("Maximum contour step size:",40), dsMax
  write(iow,'(3X,A40,I16)')adjl("Number of contour steps:",40), nsEdw
  write(6  ,'(3X,A40,I16)')adjl("Number of contour steps:",40), nsEdw
endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine ComputeContourStep
