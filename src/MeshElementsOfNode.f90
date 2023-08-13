!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine MeshElementsOfNode(numOfElementsOfNodeMax, numElementsOfNode)
!----------------------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod, only: iow
use geometry_mod,    only: numNodes, numElementsTypeDomain, elementOfNode, globalNodeIdTypeDomain
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none

integer :: element, localNodeId, globalNodeId
integer, intent(out) :: numOfElementsOfNodeMax
integer, intent(out), dimension(numNodes) :: numElementsOfNode
!----------------------------------------------------------------------------------------------------------------------------------!
numElementsOfNode = 0

numOfElementsOfNodeMax = 0
do element = 1, numElementsTypeDomain
 do localNodeId = 1, 4
   globalNodeId = globalNodeIdTypeDomain(localNodeId,element)

   numElementsOfNode(globalNodeId) = numElementsOfNode(globalNodeId) + 1

   numOfElementsOfNodeMax = MAX(numOfElementsOfNodeMax, numElementsOfNode(globalNodeId))
 enddo
enddo

write(iow,'(3X,"Maximum number of elements per node: ",I18)') numOfElementsOfNodeMax
write(6  ,'(3X,"Maximum number of elements per node: ",I18)') numOfElementsOfNodeMax

allocate(elementOfNode(numNodes, numOfElementsOfNodeMax))
elementOfNode = 0

numElementsOfNode = 0
do element = 1, numElementsTypeDomain
  do localNodeId = 1, 4
    globalNodeId = globalNodeIdTypeDomain(localNodeId,element)

    numElementsOfNode(globalNodeId) = numElementsOfNode(globalNodeId) + 1

    elementOfNode(globalNodeId, numElementsOfNode(globalNodeId)) = element
  enddo
enddo

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine MeshElementsOfNode
