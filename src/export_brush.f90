subroutine export_brush(num_gpoints, numnp, phia_gr, phia_gr_indiv, volnp, file_name, dist_from_surf)
!-----------------------------------------------------------------------------------------------------------!
use write_helper, only : adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: num_gpoints, numnp
integer             :: kk, ii

character(40), intent(in) :: file_name

real(8), intent(in), dimension(numnp)              :: phia_gr, volnp, dist_from_surf
real(8), intent(out), dimension(numnp,num_gpoints) :: phia_gr_indiv
real(8)                                            :: numer, denom, r_center_surf, msq_brush_all, &
                                                      & msq_brush_of_chain_ave, msq_brush_of_chain_std
real(8), dimension(num_gpoints)                    :: msq_brush_of_chain
!-----------------------------------------------------------------------------------------------------------!
write(6,'(6X,A40)')adjl("*mean brush thickness..",40)

msq_brush_of_chain     = 0.d0
msq_brush_all          = 0.d0
msq_brush_of_chain_ave = 0.d0
msq_brush_of_chain_std = 0.d0

open (unit=120, file = file_name)
do ii = 1, num_gpoints
    numer = 0.d0
    denom = 0.d0
    do kk = 1, numnp
        r_center_surf = dist_from_surf(kk)
        numer         = numer + r_center_surf**2 * phia_gr_indiv(kk,ii) * volnp(kk)
        denom         = denom + phia_gr_indiv(kk,ii) * volnp(kk)
    enddo
    msq_brush_of_chain(ii) = SQRT(numer / denom)
    write(120,'(I19,F19.9)') ii, msq_brush_of_chain(ii)
enddo

!compute average
do ii = 1, num_gpoints
    msq_brush_of_chain_ave = msq_brush_of_chain_ave + msq_brush_of_chain(ii)
enddo
msq_brush_of_chain_ave = msq_brush_of_chain_ave / REAL(num_gpoints)

!compute stdev
do ii = 1, num_gpoints
    msq_brush_of_chain_std = msq_brush_of_chain_std + (msq_brush_of_chain(ii) - msq_brush_of_chain_ave)**2
enddo
msq_brush_of_chain_std = SQRT(msq_brush_of_chain_std / REAL(num_gpoints))

!find the brush thickness from the total phi_gr
numer = 0.d0
denom = 0.d0
do kk = 1, numnp
    r_center_surf = dist_from_surf(kk)
    numer         = numer + r_center_surf**2 * phia_gr(kk) * volnp(kk)
    denom         = denom + phia_gr(kk) * volnp(kk)
enddo
msq_brush_all = SQRT(numer / denom)
write(120,'(A19,F19.9)') "mean",  msq_brush_of_chain_ave
write(120,'(A19,F19.9)') "stdev", msq_brush_of_chain_std
write(120,'(A19,F19.9)') "all",   msq_brush_all
close(120)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_brush
