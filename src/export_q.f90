!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_q(ns, q_final, q_type)
!--------------------------------------------------------------------!
use geometry_mod,     only: numnp, xc, ndm
use write_helper_mod, only: adjl
use constants_mod,    only: tol
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: time_step, ii, jj

character(2)  :: q_type
character(20) :: file_name

real(8), intent(in), dimension(ns+1,numnp) :: q_final
real(8)                                    :: iq_final
!--------------------------------------------------------------------!
write(6,'(2X,A23,A5,A8)')"Exporting propagator of",q_type," chains."

write(file_name,'("o.q",A2)') q_type

open(unit=363, file = file_name)
do ii = 1, numnp
    do jj = 1, ndm
        write(363,'(E20.9)',advance='no') xc(jj,ii)
    enddo
    do time_step = 1, ns+1
        iq_final = q_final(time_step,ii)
        if (DABS(iq_final)<tol) then
            iq_final = 0.d0
        endif
        write(363,'(E20.9)',advance='no') iq_final
    enddo
    write(363,*)
enddo
close(363)

return
!--------------------------------------------------------------------!
end subroutine export_q
