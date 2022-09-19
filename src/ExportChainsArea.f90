!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ExportChainsArea(nodeBelongsToDirichletFace, elemcon, cell_of_np, chain_type, rg2OfMonomer, chainLength, &
                              ns_ed, ds_ed, q_final, phi, ww)
!-----------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: segmentBulkDensity
use hist_mod,         only: nbin
use write_helper_mod, only: adjl
use constants_mod,    only: A3_to_m3
use arrays_mod,       only: nodeVolume
use geometry_mod,     only: numNodes
use delta_mod,        only: targetNumGraftedChains, graftPointId, graftPointValue
use fhash_module__ints_double
!-----------------------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------------------!
integer, intent(in), dimension(numNodes) :: cell_of_np
integer, intent(in)                      :: ns_ed
integer                                  :: bin, kk, ii, gnode_id

type(fhash_type__ints_double), intent(inout) :: elemcon

character(len=2), intent(in) :: chain_type

logical, intent(in), dimension(numNodes) :: nodeBelongsToDirichletFace
logical, dimension(numNodes)             :: nodeBelongsToDirichletFaceNew

real(8), intent(in), dimension(ns_ed+1)          :: ds_ed
real(8), intent(in), dimension(ns_ed+1,numNodes) :: q_final
real(8), intent(in), dimension(numNodes)         :: phi, ww
real(8), intent(in)                              :: rg2OfMonomer, chainLength
real(8), dimension(2,numNodes)                   :: qshape
real(8), dimension(ns_ed+1,numNodes)             :: qshape_final
real(8), dimension(nbin)                         :: p_cross, n_shape
real(8)                                          :: sum_qshape=0.d0, sum_Q=0.d0, sum_phi=0.d0
!-----------------------------------------------------------------------------------------------------------------------!
do bin = 7, 30
  write(6,*) "bin = ",bin

  call FemMatrixAssemble(rg2OfMonomer, ww)

  nodeBelongsToDirichletFaceNew = nodeBelongsToDirichletFace
  do kk = 1, numNodes
    if (cell_of_np(kk).eq.bin) nodeBelongsToDirichletFaceNew(kk) = .true.
  enddo

  qshape       = 0.0d0
  qshape_final = 0.0d0
  if (chain_type.eq."mx") then
    write(6,'(6X,A40)')adjl("Exporting matrix chains per area.",40)
    qshape(1,:)       = 1.0d0
    qshape_final(1,:) = 1.0d0
  endif
  if (chain_type.eq."gr") then
    write(6,'(6X,A40)')adjl("Exporting grafted chains per area.",40)

    do ii = 1, targetNumGraftedChains
      gnode_id = graftPointId(ii)

      qshape(1,gnode_id)       = graftPointValue(ii)
      qshape_final(1,gnode_id) = graftPointValue(ii)
    enddo
  endif

  call SolverEdwards(ds_ed, ns_ed, qshape, qshape_final, nodeBelongsToDirichletFaceNew, elemcon)

  sum_qshape = 0.0d0
  sum_Q      = 0.0d0
  sum_phi    = 0.0d0
  do kk = 1, numNodes
    sum_qshape = sum_qshape + qshape_final(ns_ed+1,kk) * nodeVolume(kk)
    sum_Q      = sum_Q      + q_final(ns_ed+1,kk)      * nodeVolume(kk)
    sum_phi    = sum_phi    + phi(kk)                  * nodeVolume(kk)
  enddo

  p_cross(bin) = 1.0d0 - sum_qshape / sum_Q

  n_shape(bin) = p_cross(bin) * segmentBulkDensity * sum_phi * A3_to_m3 / chainLength !/layer_area

  write(6,*) "p_cross: ", p_cross(bin)
  write(6,*) "n_shape: ", n_shape(bin)
enddo

return
!-----------------------------------------------------------------------------------------------------------------------!
end subroutine ExportChainsArea
