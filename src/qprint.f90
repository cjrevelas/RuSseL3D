subroutine qprint(q_final, q_type)
!--------------------------------------------------------------------!
use xdata
use constants
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: time_step, i1, ii2

character(4)                               :: q_type
character(20)                              :: file_name

real(8), intent(in), dimension(numnp,ns+1) :: q_final
real(8)                                    :: iq_final

!profile section
real(8), parameter                         :: rtol = 2.d0
integer, parameter                         :: max_time_step = 40
real(8), dimension(3)                      :: r_center
integer                                    :: d2, d3
!--------------------------------------------------------------------!
write(file_name,'(''q'',A4,''.out.txt'')') q_type

open(unit=363, file = file_name)
do i1 = 1, numnp
    do ii2 = 1, ndm
        write(363,'(E20.9)',advance='no') xc(ii2,i1)
    enddo
    do time_step = 1, ns+1
        iq_final = q_final(i1,time_step)
        if (dabs(iq_final)<tol) then
            iq_final = 0.d0
        endif
        write(363,'(E20.9)',advance='no') iq_final
    enddo
    write(363,*)
enddo
close(363)
!--------------------------------------------------------------------!
!the current section prints q's across the three perpendicular lines
!crossing the point r_center
r_center(1) = 0.d0
r_center(2) = 0.d0
r_center(3) = 0.d0

file_name = ""
write(file_name,'(''q'',A4,''_prof.out.txt'')') q_type
open(unit=363, file = file_name)

do ii2 = 1, ndm
    write(363,*) "dim ", ii2
    write(363,*) "pos times..."
    write(363,*)

    d2 = mod(ii2+0,3)+1
    d3 = mod(ii2+1,3)+1

    do i1 = 1, numnp
        if ( dabs(xc(d2,i1) - r_center(d2)) < rtol .and.  &
   &         dabs(xc(d3,i1) - r_center(d3)) < rtol ) then
            write(363,'(1(E20.9))',advance='no') xc(ii2,i1)
            !do time_step = 1, max_time_step
            do time_step = 1, ns+1
                write(363,'(E20.9)',advance='no') q_final(i1,time_step)
            enddo
            write(363,*)
        endif
    enddo
    write(363,*)
enddo

close(363)

return
!--------------------------------------------------------------------!
end subroutine qprint
