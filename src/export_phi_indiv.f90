!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_phi_indiv(num_gpoints, numnp, xc, phi_gr_indiv)
!------------------------------------------------------------------------------------------------------!
use error_handing_mod
use iofiles_mod,      only: phi_profile_indiv
use parser_vars_mod,  only: chainlen_gr, rho_mol_bulk
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: num_gpoints, numnp
integer             :: ii, kk

real(8), intent(in), dimension(3,numnp)           :: xc
real(8), intent(in), dimension(numnp,num_gpoints) :: phi_gr_indiv
real(8), dimension(num_gpoints)                   :: nch_gr
!------------------------------------------------------------------------------------------------------!
do ii = 1, num_gpoints
    call compute_number_of_chains(numnp, chainlen_gr, rho_mol_bulk, phi_gr_indiv(:,ii), nch_gr(ii))
enddo

open (unit=120, file = phi_profile_indiv)
write(120,'(4(A19))',advance='no') "nch", "#", "#", "#"
do ii = 1, num_gpoints
    write(120,'(E19.9E3)',advance='no') nch_gr(ii)
enddo

write(120,*)
write(120,'(4(A19))',advance='no') "np", "x", "y", "z"
do ii = 1, num_gpoints
    write(120,'(I19)',advance='no') ii
enddo

write(120,*)
do kk = 1, numnp
    write(120,'(I19,3(E19.9E3))',advance='no') kk, xc(1,kk), xc(2,kk), xc(3,kk)
    do ii = 1, num_gpoints
       write(120,'(E19.9E3)',advance='no') phi_gr_indiv(kk,ii)
    enddo
    write(120,*)
enddo
close(120)
!------------------------------------------------------------------------------------------------------!
end subroutine export_phi_indiv
