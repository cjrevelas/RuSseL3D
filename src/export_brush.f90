!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_brush(num_gpoints, numnp, phi_gr, phi_gr_indiv, volnp, file_name, dist_from_surf)
!-----------------------------------------------------------------------------------------------------------!
use write_helper_mod, only : adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: num_gpoints, numnp
integer             :: kk, ii

character(40), intent(in) :: file_name

real(8), intent(in), dimension(numnp)              :: phi_gr, volnp, dist_from_surf
real(8), intent(out), dimension(numnp,num_gpoints) :: phi_gr_indiv
real(8)                                            :: numer, denom, r_center_surf, msq_brush_all, &
                                                      & msq_brush_of_chain_ave, msq_brush_of_chain_std
real(8), dimension(num_gpoints)                    :: msq_brush_of_chain
!-----------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exportin mean brush thickness.",40)

msq_brush_of_chain     = 0.0d0
msq_brush_all          = 0.0d0
msq_brush_of_chain_ave = 0.0d0
msq_brush_of_chain_std = 0.0d0

! Find the brush thickness from the individual phi_gr
open (unit=120, file = file_name)
do ii = 1, num_gpoints
    numer = 0.0d0
    denom = 0.0d0

    do kk = 1, numnp
        r_center_surf = dist_from_surf(kk)
        numer         = numer + r_center_surf**2.0d0 * phi_gr_indiv(kk,ii) * volnp(kk)
        denom         = denom + phi_gr_indiv(kk,ii) * volnp(kk)
    enddo

    msq_brush_of_chain(ii) = SQRT(numer / denom)
    write(120,'(I8,2X,E16.9E3)') ii, msq_brush_of_chain(ii)
enddo

! Compute average
do ii = 1, num_gpoints
    msq_brush_of_chain_ave = msq_brush_of_chain_ave + msq_brush_of_chain(ii)
enddo
msq_brush_of_chain_ave = msq_brush_of_chain_ave / REAL(num_gpoints)

! Compute stdev
do ii = 1, num_gpoints
    msq_brush_of_chain_std = msq_brush_of_chain_std + (msq_brush_of_chain(ii) - msq_brush_of_chain_ave)**2
enddo
msq_brush_of_chain_std = SQRT(msq_brush_of_chain_std / REAL(num_gpoints))

! Find the brush thickness from the total phi_gr
numer = 0.0d0
denom = 0.0d0

do kk = 1, numnp
    r_center_surf = dist_from_surf(kk)
    numer         = numer + r_center_surf**2.0d0 * phi_gr(kk) * volnp(kk)
    denom         = denom + phi_gr(kk) * volnp(kk)
enddo
msq_brush_all = SQRT(numer / denom)

write(120,'(A8,2X,E16.9E3)') "mean",  msq_brush_of_chain_ave
write(120,'(A8,2X,E16.9E3)') "stdev", msq_brush_of_chain_std
write(120,'(A8,2X,E16.9E3)') "all",   msq_brush_all
close(120)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_brush
