!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportAdsorbed(node_belongs_to_dirichlet_face, adsorbed)
!-----------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: mumpsMatrixType, rg2OfMatrixMonomer, lengthMatrix, numEdwPointsMatrix, numConvolPointsMatrix
use write_helper_mod, only: adjl
use arrays_mod,       only: ds_mx_ed, xs_mx_ed, xs_mx_conv, coeff_mx_conv, qmx_interp_mm, phi_mx, ww
use geometry_mod,     only: numNodes, nodeCoord
!-----------------------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------------------!
integer :: ii, kk

logical, intent(in), dimension(numNodes) :: node_belongs_to_dirichlet_face, adsorbed
logical, dimension(numNodes)             :: node_belongs_to_dirichlet_face_new

real(8), dimension(2,numNodes)                       :: qfree
real(8), dimension(numEdwPointsMatrix+1,numNodes)    :: qfree_final
real(8), dimension(numConvolPointsMatrix+1,numNodes) :: qfree_interp, qads
real(8), dimension(numNodes)                         :: phi_free, phi_ads, phi_loop, phi_tail
!-----------------------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting ads vs free density profiles.",40)

call FemMatrixAssemble(rg2OfMatrixMonomer, ww)

qfree            = 0.0d0
qfree_final      = 0.0d0
qfree(1,:)       = 1.0d0
qfree_final(1,:) = 1.0d0

node_belongs_to_dirichlet_face_new = node_belongs_to_dirichlet_face

do kk = 1, numNodes
  if (adsorbed(kk)) then
    qfree(1,kk)       = 0.0d0
    qfree_final(1,kk) = 0.0d0

    node_belongs_to_dirichlet_face_new(kk) = .True.
  endif
enddo

call SolverEdwards(ds_mx_ed, numEdwPointsMatrix, mumpsMatrixType, qfree, qfree_final, node_belongs_to_dirichlet_face_new)

do ii = 1, numNodes
  call interp_linear(1, numEdwPointsMatrix+1, xs_mx_ed, qfree_final(:,ii), numConvolPointsMatrix+1, xs_mx_conv, qfree_interp(:,ii))
enddo

do ii = 1, numConvolPointsMatrix+1
  do kk = 1, numNodes
    qads(ii,kk) = qmx_interp_mm(ii,kk) - qfree_interp(ii,kk)
  enddo
enddo

! Determine the density profiles
call ContourConvolution(numNodes, lengthMatrix, numConvolPointsMatrix, coeff_mx_conv, qfree_interp, qfree_interp, phi_free)
call ContourConvolution(numNodes, lengthMatrix, numConvolPointsMatrix, coeff_mx_conv, qads, qads, phi_loop)
call ContourConvolution(numNodes, lengthMatrix, numConvolPointsMatrix, coeff_mx_conv, qfree_interp, qads, phi_tail)

open(unit=125, file="o.ads_free_mx")
write(125,'(9(2X,A16))') "kk", 'x', 'y', 'z', "phi_free", "phi_ads", "phi_loop", "phi_tail", "phi_mx"
do kk = 1, numNodes
  phi_ads(kk) = phi_mx(kk) - phi_free(kk)
  write(125,'(2X,I16,8(2X,E19.9E2))') kk, nodeCoord(1,kk), nodeCoord(2,kk), nodeCoord(3,kk), phi_free(kk), phi_ads(kk), phi_loop(kk), phi_tail(kk), phi_mx(kk)
enddo
close(125)

return
!-----------------------------------------------------------------------------------------------------------------------!
end subroutine ExportAdsorbed
