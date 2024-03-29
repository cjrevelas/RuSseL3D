subroutine ExportVtuProfiles(phi_mx, phi_gr, ww)
!----------------------------------------------------------------------------------------------------------------------------------!
use geometry_mod, only: nodeCoord, numNodes, numElementsTypeDomain, numDimensions, numNodesLocalTypeDomain, globalNodeIdTypeDomain
use iofiles_mod, only: IO_VtuProfiles
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
real(8), intent(in), dimension(numNodes) :: phi_mx, phi_gr, ww

integer, allocatable, dimension(:) :: offset
integer                            :: ii, jj
integer                            :: vtuTetType = 10

character(len=6)  :: phi_matrix  = "phi_mx"
character(len=6)  :: phi_grafted = "phi_gr"
character(len=5)  :: field       = "field"
character(len=38) :: line1       = '<?xml version="1.0" encoding="UTF-8"?>'
character(len=74) :: line2       = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
character(len=23) :: line3a      = '<Piece NumberOfPoints="'
character(len=17) :: line3b      = '" NumberOfCells="'
character(len=2)  :: line3c      = '">'
character(len=32) :: line4a      = '<DataArray type="Float64" Name="'
character(len=40) :: line4b      = '" Format="ascii" NumberOfComponents="1">'
character(len=83) :: line5       = '<DataArray type="Float64" Name="coordinates" Format="ascii" NumberOfComponents="3">'
character(len=82) :: line6       = '<DataArray type="Int32" Name="connectivity" Format="ascii" NumberOfComponents="1">'
character(len=77) :: line7       = '<DataArray type="Int32" Name="offsets" Format="ascii" NumberOfComponents="1">'
character(len=75) :: line8       = '<DataArray type="UInt8" Name="types" Format="ascii" NumberOfComponents="1">'
!----------------------------------------------------------------------------------------------------------------------------------!
if (numNodesLocalTypeDomain.eq.4) then
  continue
else
  return
endif

open(unit=1111, file=IO_VtuProfiles)
write(1111,'(A38)') line1
write(1111,'(A74)') line2
write(1111,'(A18)') "<UnstructuredGrid>"
write(1111,'(A23,1X,I6,A17,1X,I8,A2)') line3a, numNodes, line3b, numElementsTypeDomain, line3c
write(1111,*)

! Export density profile of matrix chains
write(1111,'(A11)') "<PointData>"
write(1111,'(A32,A6,A40)') line4a, phi_matrix, line4b
do ii = 1, numNodes
  write(1111,'(F19.10)') phi_mx(ii)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,*)

! Export density profile of grafted chains
write(1111,'(A32,A6,A40)') line4a, phi_grafted, line4b
do ii = 1, numNodes
  write(1111,'(F19.10)') phi_gr(ii)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,*)

! Export field configuration
write(1111,'(A32,A5,A40)') line4a, field, line4b
do ii = 1, numNodes
  write(1111,'(F19.10)') ww(ii)
enddo
write(1111,'(A12)') "</DataArray>"
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
write(1111,'(A8)') "</Piece>"
write(1111,'(A19)') "</UnstructuredGrid>"
write(1111,'(A10)') "</VTKFile>"
close(1111)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine ExportVtuProfiles
