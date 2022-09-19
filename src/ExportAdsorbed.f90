!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportAdsorbed(nodeBelongsToDirichletFace, elemcon, adsorbed)
!-----------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: rg2OfMatrixMonomer, lengthMatrix, &
                            numEdwPointsMatrix, numConvolPointsMatrix
use write_helper_mod, only: adjl
use arrays_mod,       only: ds_mx_ed, xs_mx_ed, xs_mx_conv, coeff_mx_conv, &
                            qqMatrixInterp, phiMatrix, wwField
use geometry_mod,     only: numNodes, nodeCoord
use fhash_module__ints_double
!-----------------------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------------------!
integer :: ii, kk

logical, intent(in), dimension(numNodes) :: nodeBelongsToDirichletFace, adsorbed
logical, dimension(numNodes)             :: nodeBelongsToDirichletFaceNew

type(fhash_type__ints_double), intent(inout) :: elemcon

real(8), dimension(2,numNodes)                       :: qqFree
real(8), dimension(numEdwPointsMatrix+1,numNodes)    :: qqFreeFinal
real(8), dimension(numConvolPointsMatrix+1,numNodes) :: qqFreeInterp, qqAdsorbed
real(8), dimension(numNodes)                         :: phiFree, phiAdsorbed, phiLoop, phiTail
!-----------------------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting ads vs free density profiles.",40)

call FemMatrixAssemble(rg2OfMatrixMonomer, wwField)

qqFree           = 0.0d0
qqFreeFinal      = 0.0d0
qqFree(1,:)      = 1.0d0
qqFreeFinal(1,:) = 1.0d0

nodeBelongsToDirichletFaceNew = nodeBelongsToDirichletFace

do kk = 1, numNodes
  if (adsorbed(kk)) then
    qqFree(1,kk)      = 0.0d0
    qqFreeFinal(1,kk) = 0.0d0

    nodeBelongsToDirichletFaceNew(kk) = .True.
  endif
enddo

call SolverEdwards(ds_mx_ed, numEdwPointsMatrix, qqFree, qqFreeFinal, nodeBelongsToDirichletFaceNew, elemcon)

do ii = 1, numNodes
  call interp_linear(1, numEdwPointsMatrix+1, xs_mx_ed, qqFreeFinal(:,ii), numConvolPointsMatrix+1, xs_mx_conv, qqFreeInterp(:,ii))
enddo

do ii = 1, numConvolPointsMatrix+1
  do kk = 1, numNodes
    qqAdsorbed(ii,kk) = qqMatrixInterp(ii,kk) - qqFreeInterp(ii,kk)
  enddo
enddo

! Determine the density profiles
call ContourConvolution(lengthMatrix, numConvolPointsMatrix, coeff_mx_conv, qqFreeInterp, qqFreeInterp, phiFree)
call ContourConvolution(lengthMatrix, numConvolPointsMatrix, coeff_mx_conv, qqAdsorbed, qqAdsorbed, phiLoop)
call ContourConvolution(lengthMatrix, numConvolPointsMatrix, coeff_mx_conv, qqFreeInterp, qqAdsorbed, phiTail)

open(unit=125, file="o.ads_free_mx")
write(125,'(9(2X,A16))') "kk", 'x', 'y', 'z', "phi_free", "phi_ads", "phi_loop", "phi_tail", "phi_mx"
do kk = 1, numNodes
  phiAdsorbed(kk) = phiMatrix(kk) - phiFree(kk)
  write(125,'(2X,I16,8(2X,E19.9E2))') kk, nodeCoord(1,kk), nodeCoord(2,kk), nodeCoord(3,kk), &
                                      phiFree(kk), phiAdsorbed(kk), phiLoop(kk), phiTail(kk), phiMatrix(kk)
enddo
close(125)

return
!-----------------------------------------------------------------------------------------------------------------------!
end subroutine ExportAdsorbed
