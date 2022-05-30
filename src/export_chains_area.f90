!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_chains_area(node_belongs_to_dirichlet_face, cell_of_np, chain_type, Rg2_per_mon, chainlen, &
                              ns_ed, ds_ed, q_final, phi, ww)
!-----------------------------------------------------------------------------------------------------------------------!
use parser_vars_mod,  only: mumpsMatrixType, segmentBulkDensity
use hist_mod,         only: nbin
use write_helper_mod, only: adjl
use constants_mod,    only: A3_to_m3
use arrays_mod,       only: volnp
use geometry_mod,     only: numNodes
use delta_mod,        only: targetNumGraftedChains, gp_init_value, gpid
!-----------------------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------------------!
integer, intent(in), dimension(numNodes) :: cell_of_np
integer, intent(in)                      :: ns_ed
integer                                  :: bin, kk, ii, gnode_id

character(len=2), intent(in) :: chain_type

logical, intent(in), dimension(numNodes) :: node_belongs_to_dirichlet_face
logical, dimension(numNodes)             :: node_belongs_to_dirichlet_face_new

real(8), intent(in), dimension(ns_ed+1)          :: ds_ed
real(8), intent(in), dimension(ns_ed+1,numNodes) :: q_final
real(8), intent(in), dimension(numNodes)         :: phi, ww
real(8), intent(in)                              :: Rg2_per_mon, chainlen
real(8), dimension(2,numNodes)                   :: qshape
real(8), dimension(ns_ed+1,numNodes)             :: qshape_final
real(8), dimension(nbin)                         :: p_cross, n_shape
real(8)                                          :: sum_qshape=0.d0, sum_Q=0.d0, sum_phi=0.d0
!-----------------------------------------------------------------------------------------------------------------------!
do bin = 7, 30
  write(6,*) "bin = ",bin

  call fem_matrix_assemble(Rg2_per_mon, ww)

  node_belongs_to_dirichlet_face_new = node_belongs_to_dirichlet_face
  do kk = 1, numNodes
    if (cell_of_np(kk).eq.bin) node_belongs_to_dirichlet_face_new(kk) = .true.
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
      gnode_id = gpid(ii)

      qshape(1,gnode_id)       = gp_init_value(ii)
      qshape_final(1,gnode_id) = gp_init_value(ii)
    enddo
  endif

  call solver_edwards(ds_ed, ns_ed, mumpsMatrixType, qshape, qshape_final, node_belongs_to_dirichlet_face_new)

  sum_qshape = 0.0d0
  sum_Q      = 0.0d0
  sum_phi    = 0.0d0
  do kk = 1, numNodes
    sum_qshape = sum_qshape + qshape_final(ns_ed+1,kk) * volnp(kk)
    sum_Q      = sum_Q      + q_final(ns_ed+1,kk)      * volnp(kk)
    sum_phi    = sum_phi    + phi(kk)                  * volnp(kk)
  enddo

  p_cross(bin) = 1.0d0 - sum_qshape / sum_Q

  n_shape(bin) = p_cross(bin) * segmentBulkDensity * sum_phi * A3_to_m3 / chainlen !/layer_area

  write(6,*) "p_cross: ", p_cross(bin)
  write(6,*) "n_shape: ", n_shape(bin)
enddo

return
!-----------------------------------------------------------------------------------------------------------------------!
end subroutine export_chains_area
