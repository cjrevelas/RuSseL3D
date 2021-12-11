!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_phi_smear(cell_of_np, cell_vol, numnp, file_name, phia_mx, phia_gr, volnp, lbin, nbin)
!-----------------------------------------------------------------------------------------------------------!
use write_helper_mod, only: adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in)                   :: numnp, nbin
integer                               :: kk, jj, bin
integer, intent(in), dimension(numnp) :: cell_of_np

character(40) :: file_name

real(8), intent(in)                   :: lbin
real(8), intent(in), dimension(nbin)  :: cell_vol
real(8), intent(in), dimension(numnp) :: phia_mx, phia_gr, volnp
real(8)                               :: np_mass
real(8), dimension(nbin)              :: phi_mx_smear, phi_g_smear
real(8), allocatable,dimension(:)     :: mass_layer
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting smeared density profiles.",40)

phi_mx_smear = 0.d0
phi_g_smear  = 0.d0

allocate(mass_layer(nbin))

mass_layer = 0.d0
do kk = 1, numnp
    bin             = cell_of_np(kk)
    np_mass         = phia_gr(kk) * volnp(kk)
    mass_layer(bin) = mass_layer(bin) + np_mass
enddo

do jj = 1, nbin
    if (cell_vol(jj) > 0) phi_g_smear(jj) = mass_layer(jj) / cell_vol(jj)
enddo

mass_layer = 0.d0
do kk = 1, numnp
    bin             = cell_of_np(kk)
    np_mass         = phia_mx(kk) * volnp(kk)
    mass_layer(bin) = mass_layer(bin) + np_mass
enddo

do jj = 1, nbin
    if (cell_vol(jj) > 0) phi_mx_smear(jj) = mass_layer(jj) / cell_vol(jj)
enddo

open (unit=120, file = file_name)
write(120,'(3A19)') "r","phi_m","phi_g"
do jj = 1, nbin
    write(120,'(F19.9,2(E19.9E3))')  (REAL(jj)-0.5d0)*lbin, phi_mx_smear(jj), phi_g_smear(jj)
enddo
close(120)
deallocate(mass_layer)

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_phi_smear
