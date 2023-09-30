!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportPropagator(numSteps, qqFinal, chainType)
!--------------------------------------------------------------------!
use geometry_mod,     only: numNodes, nodeCoord, numDimensions
use write_helper_mod, only: adjl
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer, intent(in) :: numSteps
integer             :: step, node, axis

character(2), intent(in) :: chainType
character(20)            :: fileName

real(8), intent(in), dimension(numSteps+1,numNodes) :: qqFinal

real(8) :: iqFinal, tol = 1.0d-12
!--------------------------------------------------------------------!
write(6,'(2X,A23,A5,A8)') "Exporting propagator of", chainType, " chains."

write(fileName,'("o.q",A2)') chainType

open(unit=363, file = fileName)
do node = 1, numNodes
  do axis = 1, numDimensions
    write(363,'(E20.9)',advance='no') nodeCoord(axis,node)
  enddo
  do step = 1, numSteps + 1
    iqFinal = qqFinal(step,node)
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
