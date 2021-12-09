subroutine export_brush99(cell_of_np, num_gpoints, numnp, file_name, phia_gr, phia_gr_indiv, volnp, lbin, nbin)
!-----------------------------------------------------------------------------------------------------------!
use write_helper, only : adjl
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer, intent(in)                   :: num_gpoints, numnp, nbin
integer, intent(in), dimension(numnp) :: cell_of_np
integer                               :: kk, ii, jj, bin

character(40) :: file_name

real(8), intent(in)                                :: lbin
real(8), intent(in), dimension(numnp)              :: phia_gr, volnp
real(8), intent(out), dimension(numnp,num_gpoints) :: phia_gr_indiv
real(8)                                            :: mass_cumulative, np_mass
real(8)                                            :: mass_total, brush99_all
real(8)                                            :: brush99_of_chain_ave, brush99_of_chain_std
real(8), dimension(num_gpoints)                    :: brush99_of_chain
real(8), allocatable, dimension(:)                 :: mass_layer
!-----------------------------------------------------------------------------------------------------------!
write(6,'(6X,A40)')adjl("*99% brush thickness..",40)

allocate(mass_layer(nbin))
open (unit=120, file = file_name)
!find the brush thickness from the total phi_gr
do jj = 1, num_gpoints
    mass_layer = 0.d0
    mass_total = 0.d0
    do kk = 1, numnp
        bin             = cell_of_np(kk)
        np_mass         = phia_gr_indiv(kk,jj) * volnp(kk)
        mass_layer(bin) = mass_layer(bin) + np_mass
        mass_total      = mass_total + np_mass
    enddo
    mass_cumulative = 0.d0
    do ii = 1, nbin
        mass_cumulative = mass_cumulative + mass_layer(ii)
        !write(*,*)jj, ii, mass_layer(ii), mass_cumulative
        if (mass_cumulative > mass_total*0.99d0) then
            brush99_of_chain(jj) = lbin * (REAL(ii) - 0.5d0)
            exit
        endif
    enddo
    write(120,'(I19,F19.9)') jj, brush99_of_chain(jj)
enddo

!compute average
do ii = 1, num_gpoints
    brush99_of_chain_ave = brush99_of_chain_ave + brush99_of_chain(ii)
enddo
brush99_of_chain_ave = brush99_of_chain_ave / REAL(num_gpoints)

!compute stdev
do ii = 1, num_gpoints
    brush99_of_chain_std = brush99_of_chain_std + (brush99_of_chain(ii) - brush99_of_chain_ave)**2
enddo
brush99_of_chain_std = SQRT(brush99_of_chain_std / REAL(num_gpoints))

!find the brush thickness from the total phi_gr
mass_layer = 0.d0
mass_total = 0.d0
do kk = 1, numnp
    bin             = cell_of_np(kk)
    np_mass         = phia_gr(kk) * volnp(kk)
    mass_layer(bin) = mass_layer(bin) + np_mass
    mass_total      = mass_total + np_mass
enddo
mass_cumulative = 0.d0
do ii = 1, nbin
    mass_cumulative = mass_cumulative + mass_layer(ii)
    !write(*,*)ii, mass_layer(ii), mass_cumulative
    if (mass_cumulative > mass_total*0.99d0) then
        brush99_all = lbin * (REAL(ii) - 0.5d0)
        exit
    endif
enddo

write(120,'(A19,F19.9)') "mean" , brush99_of_chain_ave
write(120,'(A19,F19.9)') "stdev", brush99_of_chain_std
write(120,'(A19,F19.9)') "all"  , brush99_all

close(120)
deallocate(mass_layer)
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_brush99