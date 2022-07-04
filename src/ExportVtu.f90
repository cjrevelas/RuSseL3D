subroutine ExportVtu(uu, chainType)
!----------------------------------------------------------------------------------------------------------------------------------!
use geometry_mod, only: nodeCoord, numNodes, numElementsTypeDomain, numDimensions, numNodesLocalTypeDomain, globalNodeIdTypeDomain
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
real(8), intent(in), dimension(numNodes) :: uu

integer, allocatable, dimension(:) :: offset
integer                            :: ii, jj
integer                            :: vtuTetType = 10

character(len=2), intent(in) :: chainType
character(len=20)            :: fileName
character(len=38)            :: line1  = '<?xml version="1.0" encoding="UTF-8"?>'
character(len=74)            :: line2  = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
character(len=23)            :: line3a = '<Piece NumberOfPoints="'
character(len=17)            :: line3b = '" NumberOfCells="'
character(len=2)             :: line3c = '">'
character(len=69)            :: line4  = '<DataArray type="Float64" Name="Dependent_variable_u" Format="ascii">'
character(len=64)            :: line5  = '<DataArray type="Float64" NumberOfComponents="3" Format="ascii">'
character(len=59)            :: line6  = '<DataArray type="Int32" Name="connectivity" Format="ascii">'
character(len=54)            :: line7  = '<DataArray type="Int32" Name="offsets" Format="ascii">'
character(len=52)            :: line8  = '<DataArray type="UInt8" Name="types" Format="ascii">'
!----------------------------------------------------------------------------------------------------------------------------------!
write(fileName,'("o.phi.",A2,".vtu")') chainType

if (numNodesLocalTypeDomain.eq.4) then
  continue
else
  return
endif

open(unit=1111, file=fileName)
write(1111,'(A38)') line1
write(1111,'(A74)') line2
write(1111,'(A18)') "<UnstructuredGrid>"
write(1111,'(A23,1X,I6,A17,1X,I8,A2)') line3a, numNodes, line3b, numElementsTypeDomain, line3c
write(1111,*)
write(1111,'(A11)') "<PointData>"
write(1111,'(A69)') line4
do ii = 1, numNodes
  write(1111,'(F18.16)') uu(ii)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,'(A12)') "</PointData>"
write(1111,*)
write(1111,'(A11)') "<CellData/>"
write(1111,'(A8)') "<Points>"
write(1111,'(A64)') line5
do ii = 1, numNodes
  write(1111,'(3(F19.15,1X))') (nodeCoord(jj,ii), jj = 1, numDimensions)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,'(A9)') "</Points>"
write(1111,*)
write(1111,'(A7)') "<Cells>"
write(1111,'(A59)') line6
do ii = 1, numElementsTypeDomain
  write(1111,'(4(I5,1X))') (globalNodeIdTypeDomain(jj,ii) - 1, jj = 1, numNodesLocalTypeDomain)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,*)
write(1111,'(A54)') line7
allocate(offset(4*numElementsTypeDomain))
offset = 0
jj     = 0
do ii = 1, 4*numElementsTypeDomain
  jj = jj + numNodesLocalTypeDomain
  offset(ii) = jj
enddo
do ii = 0, numElementsTypeDomain-1
  write(1111,'(4(I7,1X))') (offset(ii*4+jj), jj = 1, numNodesLocalTypeDomain)
enddo
deallocate(offset)
write(1111,'(A12)') "</DataArray>"
write(1111,*)
write(1111,'(A52)') line8
!add tet types
do ii = 0, numElementsTypeDomain-numNodesLocalTypeDomain*numNodesLocalTypeDomain, numNodesLocalTypeDomain*numNodesLocalTypeDomain
  write(1111,'(16(I2,1X))') (vtuTetType, jj = 1, numNodesLocalTypeDomain*numNodesLocalTypeDomain)
enddo
write(1111,'(9(I2,1X))') (vtuTetType, jj = 1, MOD(numElementsTypeDomain,numNodesLocalTypeDomain*numNodesLocalTypeDomain))
write(1111,'(A12)') "</DataArray>"
write(1111,'(A8)') "</Cells>"
write(1111,*)
write(1111,'(A8)') "</Piece>"
write(1111,'(A19)') "</UnstructuredGrid>"
write(1111,'(A10)') "</VTKFile>"
close(1111)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine ExportVtu
