!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ComputeIndivProfile(numNodes, elemcon, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, ww, &
                             targetNumGraftedChains, gpid, gp_init_value, phi_gr_indiv)
!------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: nodeBelongsToDirichletFace
use write_helper_mod, only: adjl
use parser_vars_mod,  only: numConvolPointsGrafted, numEdwPointsGrafted, lengthGrafted, &
                            mumpsMatrixType, rg2OfGraftedMonomer
use geometry_mod,     only: nodeBelongsToDirichletFace
use write_helper_mod, only: adjl
use parser_vars_mod,  only: numConvolPointsGrafted, numEdwPointsGrafted, lengthGrafted, mumpsMatrixType, rg2OfGraftedMonomer
use error_handing_mod
use fhash_module__ints_double
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                                    :: numNodes, targetNumGraftedChains
integer, intent(in), dimension(targetNumGraftedChains) :: gpid
integer                                                :: ii, jj

type(fhash_type__ints_double), intent(inout) :: elemcon

real(8), intent(in), dimension(targetNumGraftedChains)            :: gp_init_value
real(8), intent(in), dimension(numNodes)                          :: ww
real(8), intent(in), dimension(numConvolPointsGrafted+1,numNodes) :: qmx_interp_mg
real(8), intent(in), dimension(numEdwPointsGrafted+1)             :: ds_gr_ed, xs_gr_ed
real(8), intent(in), dimension(numConvolPointsGrafted+1)          :: xs_gr_conv, coeff_gr_conv
real(8), intent(out), dimension(numNodes,targetNumGraftedChains)  :: phi_gr_indiv
real(8), dimension(2,numNodes)                                    :: qgr
real(8), dimension(numEdwPointsGrafted+1,numNodes)                :: qgr_final
real(8), dimension(numConvolPointsGrafted+1,numNodes)             :: qgr_interp
!------------------------------------------------------------------------------------------------------!
write(6,'(2X,A43)')adjl("Computing indiv profiles of grafted chains.",43)

call FemMatrixAssemble(rg2OfGraftedMonomer, ww)

do ii = 1, targetNumGraftedChains
  qgr       = 0.0d0
  qgr_final = 0.0d0

  qgr(1,gpid(ii))       = gp_init_value(ii)
  qgr_final(1,gpid(ii)) = gp_init_value(ii)

  write(6, '(2x,A19,I7,A3)', advance='no') "Grafting point id: ", gpid(ii), " ->"
  call SolverEdwards(ds_gr_ed, numEdwPointsGrafted, mumpsMatrixType, qgr, qgr_final, nodeBelongsToDirichletFace, elemcon)

  do jj = 1, numNodes
    call interp_linear(1, numEdwPointsGrafted+1, xs_gr_ed, qgr_final(:,jj), numConvolPointsGrafted+1, xs_gr_conv, qgr_interp(:,jj))
  enddo

  call ContourConvolution(numNodes, lengthGrafted, numConvolPointsGrafted, coeff_gr_conv, qgr_interp, qmx_interp_mg, phi_gr_indiv(:,ii))
enddo

return
!------------------------------------------------------------------------------------------------------!
end subroutine ComputeIndivProfile
