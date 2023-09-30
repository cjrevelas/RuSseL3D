!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportSmearedProfile(cellId, cellVolume, numNodes, fileName, phiMatrix, phiGrafted, nodeVolume, binLength, numBins)
!-----------------------------------------------------------------------------------------------------------!
use write_helper_mod, only: adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in)                      :: numNodes, numBins
integer, intent(in), dimension(numNodes) :: cellId
integer                                  :: node, bin

character(40) :: fileName

real(8), intent(in)                      :: binLength
real(8), intent(in), dimension(numBins)  :: cellVolume
real(8), intent(in), dimension(numNodes) :: phiMatrix, phiGrafted, nodeVolume
real(8), dimension(numBins)              :: phiMatrixSmeared, phiGraftedSmeared
real(8), allocatable,dimension(:)        :: massLayer
real(8)                                  :: nodeMass
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting smeared density profiles.",40)

phiMatrixSmeared  = 0.0d0
phiGraftedSmeared = 0.0d0

allocate(massLayer(numBins))

massLayer = 0.0d0
do node = 1, numNodes
  bin            = cellId(node)
  nodeMass       = phiGrafted(node) * nodeVolume(node)
  massLayer(bin) = massLayer(bin) + nodeMass
enddo

do bin = 1, numBins
  if (cellVolume(bin) > 0) phiGraftedSmeared(bin) = massLayer(bin) / cellVolume(bin)
enddo

massLayer = 0.0d0
do node = 1, numNodes
  bin            = cellId(node)
  nodeMass       = phiMatrix(node) * nodeVolume(node)
  massLayer(bin) = massLayer(bin) + nodeMass
enddo

do bin = 1, numBins
  if (cellVolume(bin) > 0) phiMatrixSmeared(bin) = massLayer(bin) / cellVolume(bin)
enddo

open (unit=120, file = fileName)
write(120,'(3A19)') "r","phi_m","phi_g"
do bin = 1, numBins
  write(120,'(F19.9,2(E19.9E3))')  (REAL(bin)-0.5d0)*binLength, phiMatrixSmeared(bin), phiGraftedSmeared(bin)
enddo
close(120)

! Deallocate memory
deallocate(massLayer)

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine ExportSmearedProfile
