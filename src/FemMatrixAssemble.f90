!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine FemMatrixAssemble(rg2OfMonomer, ww)
!----------------------------------------------------------------------------------------------------------------!
use kcw_mod,      only: F_m
use geometry_mod, only: numNodes, numElementsTypeDomain, numDimensions, numNodesLocalTypeDomain, &
                        numBulkNodePairs, nodePairId, globalNodeIdTypeDomain, nodeCoord
use iofiles_mod,  only: IO_matrixAssembly
!----------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------!
integer, dimension(numNodesLocalTypeDomain) :: globalIndex
integer                                     :: element, globalNodeIndex, localNodeIndexOne, localNodeIndexTwo
integer                                     :: axis, pair, gaussPoint, numGaussPoints, gaussLevel

real(8), intent(in)                                       :: rg2OfMonomer
real(8), intent(in), dimension(numNodes)                  :: ww
real(8), dimension(numDimensions,numNodesLocalTypeDomain) :: coord
real(8), dimension(4,11)                                  :: shapeFunction
real(8), dimension(5,11)                                  :: gaussPointData
real(8)                                                   :: jacobDeterm
!----------------------------------------------------------------------------------------------------------------!
pair = 0

F_m%c = 0.0d0
F_m%k = 0.0d0
F_m%g = 0.0d0
F_m%w = 0.0d0

do element = 1, numElementsTypeDomain
  do localNodeIndexOne = 1, numNodesLocalTypeDomain
    globalIndex(localNodeIndexOne) = globalNodeIdTypeDomain(localNodeIndexOne,element)
    do axis = 1, numDimensions
      coord(axis,localNodeIndexOne) = nodeCoord(axis, globalIndex(localNodeIndexOne))
    enddo
  enddo

  ! Set up for gauss quadrature
  gaussLevel = 3
  CALL FemGaussPoints(gaussLevel, numGaussPoints, gaussPointData)

  do gaussPoint = 1, numGaussPoints

    pair = numNodesLocalTypeDomain * numNodesLocalTypeDomain * (element-1)

    CALL FemShapeFunctions(gaussPointData(1,gaussPoint), coord, numDimensions, numNodesLocalTypeDomain, jacobDeterm, shapeFunction)

    do localNodeIndexOne = 1, numNodesLocalTypeDomain
      do localNodeIndexTwo = 1, numNodesLocalTypeDomain
        globalNodeIndex = globalIndex(localNodeIndexTwo)

        pair = pair + 1

        F_m%c(pair) = F_m%c(pair) + shapeFunction(4,localNodeIndexOne) * shapeFunction(4,localNodeIndexTwo) * jacobDeterm * gaussPointData(5,gaussPoint)

        F_m%k(pair) = F_m%k(pair) + rg2OfMonomer * (shapeFunction(1,localNodeIndexOne)*shapeFunction(1,localNodeIndexTwo)  + &
                                                    shapeFunction(2,localNodeIndexOne)*shapeFunction(2,localNodeIndexTwo)  + &
                                                    shapeFunction(3,localNodeIndexOne)*shapeFunction(3,localNodeIndexTwo)) * jacobDeterm * gaussPointData(5,gaussPoint)

        F_m%w(pair) = F_m%w(pair) + ww(globalNodeIndex) * shapeFunction(4,localNodeIndexOne) * shapeFunction(4,localNodeIndexTwo) * jacobDeterm * gaussPointData(5,gaussPoint)
      enddo
    enddo
  enddo
enddo

! Assembly global matrix using element matrices and nodePairId hash matrix created in parser_mesh.f90
do pair = 1, numBulkNodePairs
  if (F_m%isZero(pair)) then
    ! Add up contributions of same pairs met multiple times
    F_m%k(nodePairId(pair)) = F_m%k(nodePairId(pair)) + F_m%k(pair)
    F_m%k(pair)             = 0.0d0
    F_m%c(nodePairId(pair)) = F_m%c(nodePairId(pair)) + F_m%c(pair)
    F_m%c(pair)             = 0.0d0
    F_m%w(nodePairId(pair)) = F_m%w(nodePairId(pair)) + F_m%w(pair)
    F_m%w(pair)             = 0.0d0
  endif
enddo

#ifdef DEBUG_OUTPUTS
open(unit=400, file = IO_matrixAssembly)
write(400,'(3(2X,A16))') "F_m%k","F_m%c","F_m%w"
do pair = 1, numBulkNodePairs
  write(400,'(3(2X,E16.9))') F_m%k(pair), F_m%c(pair), F_m%w(pair)
enddo
close(400)
#endif

return
!----------------------------------------------------------------------------------------------------------------!
end subroutine FemMatrixAssemble
