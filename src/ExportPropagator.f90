!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportPropagator(ns, qFinal, chainType)
!--------------------------------------------------------------------!
use geometry_mod,     only: numNodes, nodeCoord, numDimensions
use write_helper_mod, only: adjl
use constants_mod,    only: tol
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: time_step, ii, jj

character(2), intent(in) :: chainType
character(20)            :: fileName

real(8), intent(in), dimension(ns+1,numNodes) :: qFinal
real(8)                                       :: iqFinal
!--------------------------------------------------------------------!
write(6,'(2X,A23,A5,A8)') "Exporting propagator of", chainType, " chains."

write(fileName,'("o.q",A2)') chainType

open(unit=363, file = fileName)
do ii = 1, numNodes
  do jj = 1, numDimensions
    write(363,'(E20.9)',advance='no') nodeCoord(jj,ii)
  enddo
  do time_step = 1, ns+1
    iqFinal = qFinal(time_step,ii)
    if (DABS(iqFinal)<tol) then
      iqFinal = 0.0d0
    endif
    write(363,'(E20.9)',advance='no') iqFinal
  enddo
  write(363,*)
enddo
close(363)

return
!--------------------------------------------------------------------!
end subroutine ExportPropagator
