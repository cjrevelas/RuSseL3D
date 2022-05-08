subroutine export_vtu(uu)
!----------------------------------------------------------------------------------------------------------------------------------!
use geometry_mod, only: xc, numnp, numel, ndm, nel, global_node_id_type_domain
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
real(8), intent(in), dimension(numnp) :: uu

integer, allocatable, dimension(:) :: offset
integer :: ii, jj
integer :: vtu_tet_type = 10

character(len=38) :: line1  = '<?xml version="1.0" encoding="UTF-8"?>'
character(len=73) :: line2  = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
character(len=23) :: line3a = '<Piece NumberOfPoints="'
character(len=17) :: line3b = '" NumberOfCells="'
character(len=2)  :: line3c = '">'
character(len=69) :: line4  = '<DataArray type="Float64" Name="Dependent_variable_u" Format="ascii">'
character(len=64) :: line5  = '<DataArray type="Float64" NumberOfComponents="3" Format="ascii">'
character(len=59) :: line6  = '<DataArray type="Int32" Name="connectivity" Format="ascii">'
character(len=54) :: line7  = '<DataArray type="Int32" Name="offsets" Format="ascii">'
character(len=52) :: line8  = '<DataArray type="UInt8" Name="types" Format="ascii">'
!----------------------------------------------------------------------------------------------------------------------------------!
open(unit=1111, file="o.phi.vtu")
write(1111,'(A38)') line1
write(1111,'(A74)') line2
write(1111,'(A18)') "<UnstructuredGrid>"
write(1111,'(A23,1X,I6,A17,1X,I8,A2)') line3a, numnp, line3b, numel, line3c
write(1111,*)
write(1111,'(A11)') "<PointData>"
write(1111,'(A69)') line4
do ii = 1, numnp
  write(1111,'(F18.16)') uu(ii)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,'(A12)') "</PointData>"
write(1111,*)
write(1111,'(A11)') "<CellData/>"
write(1111,'(A8)') "<Points>"
write(1111,'(A64)') line5
do ii = 1, numnp
  write(1111,'(3(F19.15,1X))') (xc(jj,ii), jj = 1, ndm)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,'(A9)') "</Points>"
write(1111,*)
write(1111,'(A7)') "<Cells>"
write(1111,'(A59)') line6
do ii = 1, numel
  write(1111,'(4(I5,1X))') (global_node_id_type_domain(jj,ii) - 1, jj = 1, nel)
enddo
write(1111,'(A12)') "</DataArray>"
write(1111,*)
write(1111,'(A54)') line7
allocate(offset(numel))
offset = 0
jj     = 0
do ii = 1, numel
  jj = jj + nel
  offset(ii) = jj
enddo
do ii = 0, numel-nel, nel
  write(1111,'(4(I5,1X))') (offset(ii+jj), jj = 1, nel)
enddo
deallocate(offset)
write(1111,*)
write(1111,'(A12)') "</DataArray>"
write(1111,*)
write(1111,'(A52)') line8
close(1111)

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine export_vtu
