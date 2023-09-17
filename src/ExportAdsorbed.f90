!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportAdsorbed(nodeBelongsToDirichletFace, elemcon, adsorbed)
!-----------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: rg2OfMatrixMonomer, lengthMatrix, &
                            numEdwPointsMatrix, numConvolPointsMatrix
use write_helper_mod, only: adjl
use arrays_mod,       only: dsEdwMatrix, xsEdwMatrix, xsConvMatrix, coeffConvMatrix, &
                            qqMatrixInterp, phiMatrix, wwField
use geometry_mod,     only: numNodes, nodeCoord
use fhash_module__ints_double
!-----------------------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------------------!
integer :: step, node

logical, intent(in), dimension(numNodes) :: nodeBelongsToDirichletFace, adsorbed
logical, dimension(numNodes)             :: nodeBelongsToDirichletFaceNew

type(fhash_type__ints_double), intent(inout) :: elemcon

real(8), dimension(2,numNodes)                       :: qqFree
real(8), dimension(numEdwPointsMatrix+1,numNodes)    :: qqFreeFinal
real(8), dimension(numConvolPointsMatrix+1,numNodes) :: qqFreeInterp, qqAdsorbed
real(8), dimension(numNodes)                         :: phiFree, phiAdsorbed, phiLoop, phiTail
!-----------------------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting ads vs free density profiles.",40)

CALL FemMatrixAssemble(rg2OfMatrixMonomer, wwField)

qqFree           = 0.0d0
qqFreeFinal      = 0.0d0
qqFree(1,:)      = 1.0d0
qqFreeFinal(1,:) = 1.0d0

nodeBelongsToDirichletFaceNew = nodeBelongsToDirichletFace

do node = 1, numNodes
  if (adsorbed(node)) then
    qqFree(1,node)      = 0.0d0
    qqFreeFinal(1,node) = 0.0d0

    nodeBelongsToDirichletFaceNew(node) = .True.
  endif
enddo

CALL SolverEdwards(dsEdwMatrix, numEdwPointsMatrix, qqFree, qqFreeFinal, nodeBelongsToDirichletFaceNew, elemcon)

do node = 1, numNodes
  CALL interp_linear(1, numEdwPointsMatrix+1, xsEdwMatrix, qqFreeFinal(:,node), numConvolPointsMatrix+1, xsConvMatrix, qqFreeInterp(:,node))
enddo

do step = 1, numConvolPointsMatrix+1
  do node = 1, numNodes
    qqAdsorbed(step,node) = qqMatrixInterp(step,node) - qqFreeInterp(step,node)
  enddo
enddo

! Determine the density profiles
CALL ContourConvolution(lengthMatrix, numConvolPointsMatrix, coeffConvMatrix, qqFreeInterp, qqFreeInterp, phiFree)
CALL ContourConvolution(lengthMatrix, numConvolPointsMatrix, coeffConvMatrix, qqAdsorbed, qqAdsorbed, phiLoop)
CALL ContourConvolution(lengthMatrix, numConvolPointsMatrix, coeffConvMatrix, qqFreeInterp, qqAdsorbed, phiTail)

open(unit=125, file="o.ads_free_mx")
write(125,'(9(2X,A16))') "node", 'x', 'y', 'z', "phi_free", "phi_ads", "phi_loop", "phi_tail", "phi_mx"
do node = 1, numNodes
  phiAdsorbed(node) = phiMatrix(node) - phiFree(node)
  write(125,'(2X,I16,8(2X,E19.9E2))') node, nodeCoord(1,node), nodeCoord(2,node), nodeCoord(3,node), &
                                      phiFree(node), phiAdsorbed(node), phiLoop(node), phiTail(node), phiMatrix(node)
enddo
close(125)

return
!-----------------------------------------------------------------------------------------------------------------------!
end subroutine ExportAdsorbed
