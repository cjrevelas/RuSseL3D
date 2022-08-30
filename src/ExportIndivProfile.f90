!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportIndivProfile(numGraftedChainsToExport, numNodes, nodeCoord, phi_gr_indiv)
!------------------------------------------------------------------------------------------------------!
use error_handing_mod
use iofiles_mod,      only: IO_indivProfile
use parser_vars_mod,  only: lengthGrafted, molarBulkDensity, exportAllGraftedChains, gpIndexToExport
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numGraftedChainsToExport, numNodes
integer             :: ii, kk

real(8), intent(in), dimension(3,numNodes)                        :: nodeCoord
real(8), intent(in), dimension(numNodes,numGraftedChainsToExport) :: phi_gr_indiv
real(8), dimension(numGraftedChainsToExport)                      :: nch_gr

integer :: gpIndex
!------------------------------------------------------------------------------------------------------!
do ii = 1, numGraftedChainsToExport
  if (exportAllGraftedChains.eq.1) then
    gpIndex = ii
  else
    gpIndex = gpIndexToExport(ii)
  endif
  call ComputeNumberOfChains(numNodes, lengthGrafted, molarBulkDensity, phi_gr_indiv(:,gpIndex), nch_gr(gpIndex))
enddo

open (unit=120, file = IO_indivProfile)
write(120,'(4(A19))',advance='no') "nch", "#", "#", "#"
do ii = 1, numGraftedChainsToExport
  if (exportAllGraftedChains.eq.1) then
    gpIndex = ii
  else
    gpIndex = gpIndexToExport(ii)
  endif
  write(120,'(E19.9E3)',advance='no') nch_gr(gpIndex)
enddo

write(120,*)
write(120,'(4(A19))',advance='no') "np", "x", "y", "z"
do ii = 1, numGraftedChainsToExport
  if (exportAllGraftedChains.eq.1) then
    gpIndex = ii
  else
    gpIndex = gpIndexToExport(ii)
  endif

  write(120,'(I19)',advance='no') gpIndex
enddo

write(120,*)
do kk = 1, numNodes
  write(120,'(I19,3(E19.9E3))',advance='no') kk, nodeCoord(1,kk), nodeCoord(2,kk), nodeCoord(3,kk)
  do ii = 1, numGraftedChainsToExport
    if (exportAllGraftedChains.eq.1) then
      gpIndex = ii
    else
      gpIndex = gpIndexToExport(ii)
    endif

    write(120,'(E19.9E3)',advance='no') phi_gr_indiv(kk,gpIndex)
  enddo
  write(120,*)
enddo
close(120)

return
!------------------------------------------------------------------------------------------------------!
end subroutine ExportIndivProfile
