!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_delta_numer(numNodes, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, &
                               coeff_gr_conv, ww, targetNumGraftedChains, gpid, delta_numer, volnp)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: nodeBelongsToDirichletFace
use parser_vars_mod,  only: numConvolPointsGrafted, numEdwPointsGrafted, lengthGrafted, &
                            mumpsMatrixType, rg2OfGraftedMonomer, molarBulkDensity
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: nodeBelongsToDirichletFace
use parser_vars_mod,  only: numConvolPointsGrafted, numEdwPointsGrafted, lengthGrafted, mumpsMatrixType, rg2OfGraftedMonomer, molarBulkDensity
use constants_mod,    only: n_avog, m3_to_A3
use write_helper_mod, only: adjl
use error_handing_mod
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                                    :: numNodes, targetNumGraftedChains
integer, intent(in), dimension(targetNumGraftedChains) :: gpid
integer                                                :: ii, jj

real(8), intent(in), dimension(numNodes)                          :: ww, volnp
real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qmx_interp_mg
real(8), intent(in), dimension(numEdwPointsGrafted+1)             :: ds_gr_ed, xs_gr_ed
real(8), intent(in), dimension(numConvolPointsGrafted+1)          :: xs_gr_conv, coeff_gr_conv
real(8), intent(out), dimension(targetNumGraftedChains)           :: delta_numer
real(8), dimension(targetNumGraftedChains)                        :: delta_anal
real(8), dimension(numNodes)                                      :: phi_gr
real(8), dimension(2,numNodes)                                    :: qgr
real(8), dimension(numEdwPointsGrafted+1,numNodes)                :: qgr_final
real(8), dimension(numConvolPointsGrafted+1,numNodes)             :: qgr_interp
real(8)                                                           :: nch_gr = 0.0d0, initValue = 0.0d0
!------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("****************************************",40)
write(6,'(2X,A40)')adjl("Updating delta of grafted chains:",40)

delta_anal = 0.0d0

! Analytic delta calculation
do ii = 1, targetNumGraftedChains
  delta_anal(ii) = 1.0d0 / volnp(gpid(ii)) * m3_to_A3
enddo

! Numerical delta calculation (i.e., through solution of the Edwards equation)
call fem_matrix_assemble(rg2OfGraftedMonomer, ww)

do ii = 1, targetNumGraftedChains
  qgr       = 0.0d0
  qgr_final = 0.0d0

  initValue = delta_anal(ii) * lengthGrafted * 1.0d0 / (qmx_interp_mg(numConvolPointsGrafted+1,gpid(ii)) * (molarBulkDensity * n_avog))

  qgr(1,gpid(ii))       = initValue
  qgr_final(1,gpid(ii)) = initValue

  write(6, '(2X,A19,I7,A3)', advance='no') "Grafting point id: ", gpid(ii), " ->"
  call solver_edwards(ds_gr_ed, numEdwPointsGrafted, mumpsMatrixType, qgr, qgr_final, nodeBelongsToDirichletFace)

  do jj = 1, numNodes
    call interp_linear(1, numEdwPointsGrafted+1, xs_gr_ed, qgr_final(:,jj), numConvolPointsGrafted+1, xs_gr_conv, qgr_interp(:,jj))
  enddo

  call contour_convolution(numNodes, lengthGrafted, numConvolPointsGrafted, coeff_gr_conv, qgr_interp, qmx_interp_mg, phi_gr)

  call compute_number_of_chains(numNodes, lengthGrafted, molarBulkDensity, phi_gr, nch_gr)

  delta_numer(ii) = delta_anal(ii) / nch_gr
enddo
write(6,'(2X,A40)')adjl("****************************************",40)
!------------------------------------------------------------------------------------------------------!
end subroutine compute_delta_numer
