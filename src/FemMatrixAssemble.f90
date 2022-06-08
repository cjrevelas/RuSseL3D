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
integer                                     :: lint, elem
integer                                     :: ii, jj, kk, ll, mm, nn, pp

real(8), intent(in)                                       :: rg2OfMonomer
real(8), intent(in), dimension(numNodes)                  :: ww
real(8), dimension(numDimensions,numNodesLocalTypeDomain) :: xl
real(8), dimension(4,11)                                  :: shp
real(8), dimension(5,11)                                  :: sv
real(8)                                                   :: xsj
!----------------------------------------------------------------------------------------------------------------!
kk = 0

F_m%c = 0.0d0
F_m%k = 0.0d0
F_m%g = 0.0d0
F_m%w = 0.0d0

do elem = 1, numElementsTypeDomain
  do ii = 1, numNodesLocalTypeDomain
    globalIndex(ii) = globalNodeIdTypeDomain(ii,elem)
    do jj = 1, numDimensions
      xl(jj,ii) = nodeCoord(jj, globalIndex(ii))
    enddo
  enddo

  ! Set up for gauss quadrature
  ll = 3
  call FemGaussPoints(ll, lint, sv)

  do ll = 1, lint

    kk = numNodesLocalTypeDomain * numNodesLocalTypeDomain * (elem-1)

    call FemShapeFunctions(sv(1,ll), xl, numDimensions, numNodesLocalTypeDomain, xsj, shp)

    do mm = 1, numNodesLocalTypeDomain
      do nn = 1, numNodesLocalTypeDomain
        pp = globalIndex(nn)

        kk = kk + 1

        F_m%c(kk) = F_m%c(kk) + shp(4,nn)*shp(4,mm)*xsj*sv(5,ll)

        F_m%k(kk) = F_m%k(kk) + rg2OfMonomer * (shp(1,nn)*shp(1,mm)+shp(2,nn)*shp(2,mm)+shp(3,nn)*shp(3,mm))*xsj*sv(5,ll)

        F_m%w(kk) = F_m%w(kk) + ww(pp)*shp(4,nn)*shp(4,mm)*xsj*sv(5,ll)
      enddo
    enddo
  enddo
enddo

! Assembly global matrix using element matrices and nodePairId hash matrix created in parser_mesh.f90
do kk = 1, numBulkNodePairs
  if (F_m%is_zero(kk)) then
    ! Add up contributions of same pairs met multiple times
    F_m%k(nodePairId(kk)) = F_m%k(nodePairId(kk)) + F_m%k(kk)
    F_m%k(kk)             = 0.0d0
    F_m%c(nodePairId(kk)) = F_m%c(nodePairId(kk)) + F_m%c(kk)
    F_m%c(kk)             = 0.0d0
    F_m%w(nodePairId(kk)) = F_m%w(nodePairId(kk)) + F_m%w(kk)
    F_m%w(kk)             = 0.0d0
  endif
enddo

#ifdef DEBUG_OUTPUTS
open(unit=400, file = IO_matrixAssembly)
write(400,'(3(2X,A16))') "F_m%k","F_m%c","F_m%w"
do kk = 1, numBulkNodePairs
  write(400,'(3(2X,E16.9))') F_m%k(kk), F_m%c(kk), F_m%w(kk)
enddo
close(400)
#endif

return
!----------------------------------------------------------------------------------------------------------------!
end subroutine FemMatrixAssemble
