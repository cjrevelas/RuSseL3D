subroutine ExportVtuIndiv(numGraftedChainsToExport, phi_gr_indiv)
!----------------------------------------------------------------------------------------------------------------------------------!
use geometry_mod,    only: nodeCoord, numNodes, numElementsTypeDomain, numDimensions, numNodesLocalTypeDomain, globalNodeIdTypeDomain
use delta_mod,       only: targetNumGraftedChains, graftPointId
use iofiles_mod,     only: IO_VtuGrafted
use parser_vars_mod, only: exportAllGraftedChains, gpIndexToExport
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
integer, intent(in) :: numGraftedChainsToExport

real(8), intent(in), dimension(numNodes, targetNumGraftedChains) :: phi_gr_indiv

integer, allocatable, dimension(:) :: offset
integer                            :: ii, jj, gpIndex, gpCounter, gpId
integer                            :: vtuTetType = 10

logical :: exportAll

character(len=25) :: chainId
character(len=38) :: line1  = '<?xml version="1.0" encoding="UTF-8"?>'
character(len=74) :: line2  = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
character(len=23) :: line3a = '<Piece NumberOfPoints="'
character(len=17) :: line3b = '" NumberOfCells="'
character(len=2)  :: line3c = '">'
character(len=32) :: line4a = '<DataArray type="Float64" Name="'
character(len=40) :: line4b = '" Format="ascii" NumberOfComponents="1">'
character(len=83) :: line5  = '<DataArray type="Float64" Name="coordinates" Format="ascii" NumberOfComponents="3">'
character(len=82) :: line6  = '<DataArray type="Int32" Name="connectivity" Format="ascii" NumberOfComponents="1">'
character(len=77) :: line7  = '<DataArray type="Int32" Name="offsets" Format="ascii" NumberOfComponents="1">'
character(len=75) :: line8  = '<DataArray type="UInt8" Name="types" Format="ascii" NumberOfComponents="1">'
!----------------------------------------------------------------------------------------------------------------------------------!
if (numNodesLocalTypeDomain.eq.4) then
  continue
else
  return
endif

open(unit=1111, file=IO_VtuGrafted)
write(1111,'(A38)') line1
write(1111,'(A74)') line2
write(1111,'(A18)') "<UnstructuredGrid>"
write(1111,'(A23,1X,I6,A17,1X,I8,A2)') line3a, numNodes, line3b, numElementsTypeDomain, line3c
write(1111,*)

if (exportAllGraftedChains.eq.1) then
  gpCounter = targetNumGraftedChains
  exportAll = .True.
elseif (exportAllGraftedChains.eq.0) then
  gpCounter = numGraftedChainsToExport
  exportAll = .False.
else
  exportAll = .False.
endif

write(1111,'(A11)') "<PointData>"
if ((exportAllGraftedChains.eq.0).or.(exportAllGraftedChains.eq.1)) then
  do ii = 1, gpCounter
    ! Export density profile of current grafted chain
    if (exportAll) then
      gpIndex = ii
    else
      gpIndex = gpIndexToExport(ii)
    endif

    gpId = graftPointId(gpIndex)

    write(chainId,'("index:",1X,I4,1X,", id:",1X,I7)') gpIndex, gpId

    write(1111,'(A32,A25,A40)') line4a, chainId, line4b
    do jj = 1, numNodes
      write(1111,'(F19.10)') phi_gr_indiv(jj,gpIndex)
    enddo
    write(1111,'(A12)') "</DataArray>"
    write(1111,*)
  enddo
elseif (exportAllGraftedChains.eq.2) then
  write(chainId,'("index:",1X,I4,1X,", id:",1X,I7)') 1, 1

  write(1111,'(A32,A25,A40)') line4a, chainId, line4b
  do jj = 1, numNodes
    write(1111,'(F19.10)') phi_gr_indiv(jj,1)
  enddo
  write(1111,'(A12)') "</DataArray>"
  write(1111,*)
endif

write(1111,'(A12)') "</PointData>"
write(1111,*)

! Export node coordinates
write(1111,'(A8)') "<Points>"
write(1111,'(A83)') line5
do ii = 1, numNodes
  write(1111,'(3(F19.10,1X))') (nodeCoord(jj,ii), jj = 1, numDimensions)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,'(A9)') "</Points>"
write(1111,*)

! Export node connectivity
write(1111,'(A7)') "<Cells>"
write(1111,'(A82)') line6
do ii = 1, numElementsTypeDomain
  write(1111,'(4(I8,1X))') (globalNodeIdTypeDomain(jj,ii) - 1, jj = 1, numNodesLocalTypeDomain)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,*)

! Export element offset
write(1111,'(A77)') line7
allocate(offset(4*numElementsTypeDomain))
offset = 0
jj     = 0
do ii = 1, 4*numElementsTypeDomain
  jj = jj + numNodesLocalTypeDomain
  offset(ii) = jj
enddo
do ii = 0, numElementsTypeDomain-1
  write(1111,'(4(I8,1X))') (offset(ii*4+jj), jj = 1, numNodesLocalTypeDomain)
enddo
deallocate(offset)
write(1111,'(A12)') "</DataArray>"
write(1111,*)

! Export element types
write(1111,'(A75)') line8
do ii = 0, numElementsTypeDomain-numNodesLocalTypeDomain*numNodesLocalTypeDomain, numNodesLocalTypeDomain*numNodesLocalTypeDomain
  write(1111,'(16(I2,1X))') (vtuTetType, jj = 1, numNodesLocalTypeDomain*numNodesLocalTypeDomain)
enddo
write(1111,'(9(I2,1X))') (vtuTetType, jj = 1, MOD(numElementsTypeDomain,numNodesLocalTypeDomain*numNodesLocalTypeDomain))
write(1111,'(A12)') "</DataArray>"
write(1111,'(A8)') "</Cells>"
write(1111,*)

! Close xml file
write(1111,'(A8)')  "</Piece>"
write(1111,'(A19)') "</UnstructuredGrid>"
write(1111,'(A10)') "</VTKFile>"
close(1111)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine ExportVtuIndiv
