!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportPropagator(ns, q_final, q_type)
!--------------------------------------------------------------------!
use geometry_mod,     only: numNodes, nodeCoord, numDimensions
use write_helper_mod, only: adjl
use constants_mod,    only: tol
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: time_step, ii, jj

character(2)  :: q_type
character(20) :: file_name

real(8), intent(in), dimension(ns+1,numNodes) :: q_final
real(8)                                       :: iq_final
!--------------------------------------------------------------------!
write(6,'(2X,A23,A5,A8)') "Exporting propagator of", q_type, " chains."

write(file_name,'("o.q",A2)') q_type

open(unit=363, file = file_name)
do ii = 1, numNodes
  do jj = 1, numDimensions
    write(363,'(E20.9)',advance='no') nodeCoord(jj,ii)
  enddo
  do time_step = 1, ns+1
    iq_final = q_final(time_step,ii)
    if (DABS(iq_final)<tol) then
      iq_final = 0.0d0
    endif
    write(363,'(E20.9)',advance='no') iq_final
  enddo
  write(363,*)
enddo
close(363)

return
!--------------------------------------------------------------------!
end subroutine ExportPropagator
