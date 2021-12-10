!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine get_contour_step(ds_ave, xs_crit, chainlen, ns_ed)
!------------------------------------------------------------------------------------------------------!
use constants,    only: pi
use parser_vars,  only: iow
use write_helper, only: adjl
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(out) :: ns_ed
integer              :: ns_part1, ns_part2

real(8), intent(in)  :: ds_ave, xs_crit, chainlen
real(8)              :: ds_max
!------------------------------------------------------------------------------------------------------!
if (chainlen.ge.xs_crit) then
    ns_part1 = 2 * NINT(0.5d0 * xs_crit / ds_ave)
    ds_max   = (xs_crit * (1.d0 - DCOS(pi * (DBLE(ns_part1+1)-1.d0) / (DBLE(ns_part1) * 2.d0)))) &
&            - (xs_crit * (1.d0 - DCOS(pi * (DBLE(ns_part1)  -1.d0) / (DBLE(ns_part1) * 2.d0))))
    ns_part2 = 2 * NINT(0.5d0*(chainlen - xs_crit)/ds_max)
    ns_ed    = ns_part1 + ns_part2 + 1

    write(iow,'(3X,A40)')adjl("Critical point smaller than chainlength.",40)
    write(6  ,'(3X,A40)')adjl("Critical point smaller than chainlength.",40)
    write(iow,'(3X,A40,F16.1)')adjl("Maximum contour step size:",40), ds_max
    write(6  ,'(3X,A40,F16.1)')adjl("Maximum contour step size:",40), ds_max
    write(iow,'(3X,A40,I16)')adjl("Contour steps before critical point:",40), ns_part1
    write(6  ,'(3X,A40,I16)')adjl("Contour steps before critical point:",40), ns_part1
    write(iow,'(3X,A40,I16)')adjl("Contour steps after critical point:",40), ns_part2
    write(6  ,'(3X,A40,I16)')adjl("Contour steps after critical point:",40), ns_part2
    write(iow,'(3X,A40,I16)')adjl("Total number of contour steps:",40), ns_ed
    write(6  ,'(3X,A40,I16)')adjl("Total number of contour steps:",40), ns_ed
else
    ns_ed  = CEILING( 1.d0 + 2.d0 * xs_crit / (ds_ave * pi) * ACOS(1.d0 - chainlen / xs_crit) )
    ds_max = (xs_crit * (1.d0 - DCOS(pi * (DBLE(ns_ed+1)-1.d0) / (DBLE(ns_part1) * 2.d0)))) &
&          - (xs_crit * (1.d0 - DCOS(pi * (DBLE(ns_ed)  -1.d0) / (DBLE(ns_part1) * 2.d0))))
    write(iow,'(3X,A40)')adjl("Critical point larger than chainlength.",40)
    write(6  ,'(3X,A40)')adjl("Critical point larger than chainlength.",40)
    write(iow,'(3X,A40,I16)')adjl("Number of contour steps:",40), ns_ed
    write(6  ,'(3X,A40,I16)')adjl("Number of contour steps:",40), ns_ed
endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine get_contour_step
