!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_brush99(cell_of_np, targetNumGraftedChains, numNodes, file_name, phi_gr, &
                          phi_gr_indiv, volnp, lbin, nbin)
!-----------------------------------------------------------------------------------------------------------!
use write_helper_mod, only: adjl
use parser_vars_mod,  only: iow, numNanopFaces, nanopRadiusEff
use geometry_mod,     only: boxLength
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in)                      :: targetNumGraftedChains, numNodes, nbin
integer, intent(in), dimension(numNodes) :: cell_of_np
integer                                  :: kk, ii, jj, bin

character(40) :: file_name

real(8), intent(in)                                              :: lbin
real(8), intent(in), dimension(numNodes)                         :: phi_gr, volnp
real(8), intent(out), dimension(numNodes,targetNumGraftedChains) :: phi_gr_indiv
real(8), allocatable, dimension(:)                               :: mass_layer
real(8), dimension(targetNumGraftedChains)                       :: brush99_of_chain
real(8)                                                          :: mass_cumulative, np_mass
real(8)                                                          :: mass_total, brush99_all, boxSize
real(8)                                                          :: brush99_of_chain_ave, brush99_of_chain_std
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting 99% brush thickness.",40)

brush99_of_chain_ave = 0.0d0
brush99_of_chain_std = 0.0d0

allocate(mass_layer(nbin))

open (unit=120, file = file_name)
! Find the brush thickness from the individual phi_gr
do jj = 1, targetNumGraftedChains
  mass_layer = 0.0d0
  mass_total = 0.0d0

  do kk = 1, numNodes
    bin             = cell_of_np(kk)
    np_mass         = phi_gr_indiv(kk,jj) * volnp(kk)
    mass_layer(bin) = mass_layer(bin) + np_mass
    mass_total      = mass_total + np_mass
  enddo

  mass_cumulative = 0.0d0
  do ii = 1, nbin
    mass_cumulative = mass_cumulative + mass_layer(ii)

    if (mass_cumulative > mass_total*0.99d0) then
      brush99_of_chain(jj) = lbin * (REAL(ii) - 0.5d0)
      exit
    endif
  enddo

  write(120,'(I8,2X,E16.9E3)') jj, brush99_of_chain(jj)
enddo

! Compute average
do ii = 1, targetNumGraftedChains
  brush99_of_chain_ave = brush99_of_chain_ave + brush99_of_chain(ii)
enddo
brush99_of_chain_ave = brush99_of_chain_ave / REAL(targetNumGraftedChains)

! Compute stdev
do ii = 1, targetNumGraftedChains
  brush99_of_chain_std = brush99_of_chain_std + (brush99_of_chain(ii) - brush99_of_chain_ave)**2
enddo
brush99_of_chain_std = SQRT(brush99_of_chain_std / REAL(targetNumGraftedChains))

! Find the brush thickness from the total phi_gr
mass_layer = 0.0d0
mass_total = 0.0d0
do kk = 1, numNodes
  bin             = cell_of_np(kk)
  np_mass         = phi_gr(kk) * volnp(kk)
  mass_layer(bin) = mass_layer(bin) + np_mass
  mass_total      = mass_total + np_mass
enddo

mass_cumulative = 0.0d0
do ii = 1, nbin
  mass_cumulative = mass_cumulative + mass_layer(ii)
  if (mass_cumulative > mass_total*0.99d0) then
    brush99_all = lbin * (REAL(ii) - 0.5d0)
    exit
  endif
enddo

write(120,'(A8,2X,E16.9E3)') "mean" , brush99_of_chain_ave
write(120,'(A8,2X,E16.9E3)') "stdev", brush99_of_chain_std
write(120,'(A8,2X,E16.9E3)') "all"  , brush99_all

! TODO: the following warning needs revision in the case of multiple nanoparticles
if (numNanopFaces/=0) then
  do ii = 1, numNanopFaces
    boxSize = min(boxLength(1), boxLength(2), boxLength(3))/2.d0 - nanopRadiusEff(ii)
  enddo

  if (brush99_all > (0.95d0 * boxSize)) then
    write(iow,'(2X,A52,2(A13,E11.3E3))') "WARNING: Grafted chains reach the box boundaries! ->", "Brush_h99%:", brush99_all, ", Box size:", boxSize
    write(6  ,'(2X,A52,2(A13,E11.3E3))') "WARNING: Grafted chains reach the box boundaries! ->", "Brush_h99%:", brush99_all, ", Box size:", boxSize
  endif
endif
write(6,'(2X,A40)')adjl("****************************************",40)
close(120)

! Deallocate memory
deallocate(mass_layer)

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_brush99
