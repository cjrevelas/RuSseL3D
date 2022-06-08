!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshElementsOfNode(numOfElementsOfNodeMax, numElementsOfNode)
!----------------------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod, only: iow
use geometry_mod, only: numNodes, numElementsTypeDomain, elementOfNode, globalNodeIdTypeDomain
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none

integer :: ii, jj, node
integer, intent(out) :: numOfElementsOfNodeMax
integer, intent(out), dimension(numNodes) :: numElementsOfNode
!----------------------------------------------------------------------------------------------------------------------------------!
numElementsOfNode = 0

numOfElementsOfNodeMax = 0
do ii = 1, numElementsTypeDomain
 do jj = 1, 4
   node = globalNodeIdTypeDomain(jj, ii)

   numElementsOfNode(node) = numElementsOfNode(node) + 1
   numOfElementsOfNodeMax = MAX(numOfElementsOfNodeMax, numElementsOfNode(node))
 enddo
enddo

write(iow,'(3X,"Maximum number of elements per node: ",I18)') numOfElementsOfNodeMax
write(6  ,'(3X,"Maximum number of elements per node: ",I18)') numOfElementsOfNodeMax

allocate(elementOfNode(numNodes, numOfElementsOfNodeMax))
elementOfNode = 0

numElementsOfNode = 0
do ii = 1, numElementsTypeDomain
  do jj = 1, 4
    node = globalNodeIdTypeDomain(jj,ii)

    numElementsOfNode(node) = numElementsOfNode(node) + 1
    elementOfNode(node, numElementsOfNode(node)) = ii
  enddo
enddo

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshElementsOfNode
