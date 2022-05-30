!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_phi_smeared(cell_of_np, cell_vol, numNodes, file_name, phi_mx, phi_gr, volnp, lbin, nbin)
!-----------------------------------------------------------------------------------------------------------!
use write_helper_mod, only: adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in)                      :: numNodes, nbin
integer, intent(in), dimension(numNodes) :: cell_of_np
integer                                  :: kk, jj, bin

character(40) :: file_name

real(8), intent(in)                      :: lbin
real(8), intent(in), dimension(nbin)     :: cell_vol
real(8), intent(in), dimension(numNodes) :: phi_mx, phi_gr, volnp
real(8), dimension(nbin)                 :: phi_mx_smeared, phi_gr_smeared
real(8), allocatable,dimension(:)        :: mass_layer
real(8)                                  :: np_mass
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting smeared density profiles.",40)

phi_mx_smeared = 0.0d0
phi_gr_smeared = 0.0d0

allocate(mass_layer(nbin))

mass_layer = 0.0d0
do kk = 1, numNodes
  bin             = cell_of_np(kk)
  np_mass         = phi_gr(kk) * volnp(kk)
  mass_layer(bin) = mass_layer(bin) + np_mass
enddo

do jj = 1, nbin
  if (cell_vol(jj) > 0) phi_gr_smeared(jj) = mass_layer(jj) / cell_vol(jj)
enddo

mass_layer = 0.0d0
do kk = 1, numNodes
  bin             = cell_of_np(kk)
  np_mass         = phi_mx(kk) * volnp(kk)
  mass_layer(bin) = mass_layer(bin) + np_mass
enddo

do jj = 1, nbin
  if (cell_vol(jj) > 0) phi_mx_smeared(jj) = mass_layer(jj) / cell_vol(jj)
enddo

open (unit=120, file = file_name)
write(120,'(3A19)') "r","phi_m","phi_g"
do jj = 1, nbin
  write(120,'(F19.9,2(E19.9E3))')  (REAL(jj)-0.5d0)*lbin, phi_mx_smeared(jj), phi_gr_smeared(jj)
enddo
close(120)

! Deallocate memory
deallocate(mass_layer)

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_phi_smeared
