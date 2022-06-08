!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportIndivProfile(targetNumGraftedChains, numNodes, nodeCoord, phi_gr_indiv)
!------------------------------------------------------------------------------------------------------!
use error_handing_mod
use iofiles_mod,      only: IO_indivProfile
use parser_vars_mod,  only: lengthGrafted, molarBulkDensity
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: targetNumGraftedChains, numNodes
integer             :: ii, kk

real(8), intent(in), dimension(3,numNodes)                      :: nodeCoord
real(8), intent(in), dimension(numNodes,targetNumGraftedChains) :: phi_gr_indiv
real(8), dimension(targetNumGraftedChains)                      :: nch_gr
!------------------------------------------------------------------------------------------------------!
do ii = 1, targetNumGraftedChains
  call ComputeNumberOfChains(numNodes, lengthGrafted, molarBulkDensity, phi_gr_indiv(:,ii), nch_gr(ii))
enddo

open (unit=120, file = IO_indivProfile)
write(120,'(4(A19))',advance='no') "nch", "#", "#", "#"
do ii = 1, targetNumGraftedChains
  write(120,'(E19.9E3)',advance='no') nch_gr(ii)
enddo

write(120,*)
write(120,'(4(A19))',advance='no') "np", "x", "y", "z"
do ii = 1, targetNumGraftedChains
  write(120,'(I19)',advance='no') ii
enddo

write(120,*)
do kk = 1, numNodes
  write(120,'(I19,3(E19.9E3))',advance='no') kk, nodeCoord(1,kk), nodeCoord(2,kk), nodeCoord(3,kk)
  do ii = 1, targetNumGraftedChains
    write(120,'(E19.9E3)',advance='no') phi_gr_indiv(kk,ii)
  enddo
  write(120,*)
enddo
close(120)
!------------------------------------------------------------------------------------------------------!
end subroutine ExportIndivProfile
