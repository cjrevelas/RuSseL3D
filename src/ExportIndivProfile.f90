!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportIndivProfile(numGraftedChainsToExport, numNodes, nodeCoord, phiGraftedIndiv)
!------------------------------------------------------------------------------------------------------!
use error_handling_mod
use iofiles_mod,     only: IO_indivProfile
use parser_vars_mod, only: lengthGrafted, exportAllGraftedChains, graftingPointIndexToExport
use delta_mod,       only: targetNumGraftedChains
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numGraftedChainsToExport, numNodes
integer             :: node, graftedChain, graftingPointIndex

real(8), intent(in), dimension(3,numNodes)                      :: nodeCoord
real(8), intent(in), dimension(numNodes,targetNumGraftedChains) :: phiGraftedIndiv
real(8), dimension(targetNumGraftedChains)                      :: numChainsGrafted
!------------------------------------------------------------------------------------------------------!
do graftedChain = 1, numGraftedChainsToExport
  if (exportAllGraftedChains.eq.1) then
    graftingPointIndex = graftedChain
  else
    graftingPointIndex = graftingPointIndexToExport(graftedChain)
  endif
  call ComputeNumberOfChains(lengthGrafted, phiGraftedIndiv(:,graftingPointIndex), numChainsGrafted(graftingPointIndex))
enddo

open (unit=120, file = IO_indivProfile)
write(120,'(4(A19))',advance='no') "nch", "#", "#", "#"
do graftedChain = 1, numGraftedChainsToExport
  if (exportAllGraftedChains.eq.1) then
    graftingPointIndex = graftedChain
  else
    graftingPointIndex = graftingPointIndexToExport(graftedChain)
  endif
  write(120,'(E19.9E3)',advance='no') numChainsGrafted(graftingPointIndex)
enddo

write(120,*)
write(120,'(4(A19))',advance='no') "np", "x", "y", "z"
do graftedChain = 1, numGraftedChainsToExport
  if (exportAllGraftedChains.eq.1) then
    graftingPointIndex = graftedChain
  else
    graftingPointIndex = graftingPointIndexToExport(graftedChain)
  endif

  write(120,'(I19)',advance='no') graftingPointIndex
enddo

write(120,*)
do node = 1, numNodes
  write(120,'(I19,3(E19.9E3))',advance='no') node, nodeCoord(1,node), nodeCoord(2,node), nodeCoord(3,node)
  do graftedChain = 1, numGraftedChainsToExport
    if (exportAllGraftedChains.eq.1) then
      graftingPointIndex = graftedChain
    else
      graftingPointIndex = graftingPointIndexToExport(graftedChain)
    endif

    write(120,'(E19.9E3)',advance='no') phiGraftedIndiv(node,graftingPointIndex)
  enddo
  write(120,*)
enddo
close(120)

return
!------------------------------------------------------------------------------------------------------!
end subroutine ExportIndivProfile
