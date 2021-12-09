subroutine init_chain_contour(discr, chainlen, ns, ds_ave, ds, xs, coeff)
!------------------------------------------------------------------------------------------------------!
use constants,   only: pi
use parser_vars, only: contour_discr_scheme
use flags,       only: uniform, nonuniform, symm, asymm, hybrid
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: discr, ns
integer             :: k1, ns_part1

real(8), intent(in)                   :: chainlen, ds_ave
real(8), intent(out), dimension(ns+1) :: ds, xs, coeff
real(8)                               :: xs_crit = 4.d01
!------------------------------------------------------------------------------------------------------!
if (contour_discr_scheme.eq.uniform) then
    ds = ds_ave
    do k1 = 2, ns+1
        xs(k1) = xs(k1-1) + ds_ave
    enddo
elseif (contour_discr_scheme.eq.nonuniform) then
    ds(1)=0.d0
    if (discr.eq.symm) then
        do k1 = 2, ns+1
            xs(k1) = chainlen * 0.5d0 * (1.d0 - DCOS(pi * (DBLE(k1)-1.d0) /  DBLE(ns)))
            ds(k1) = xs(k1) - xs(k1-1)
        enddo
    else if (discr.eq.asymm) then
        do k1 = 2, ns+1
            xs(k1) = chainlen *         (1.d0 - DCOS(pi * (DBLE(k1)-1.d0) / (DBLE(ns) * 2.d0)))
            ds(k1) = xs(k1) - xs(k1-1)
        enddo
    elseif (discr.eq.hybrid) then
        ns_part1 = 2 * NINT(0.5d0 * xs_crit / ds_ave)
        do k1 = 2, ns_part1+1 
            xs(k1) = xs_crit *         (1.d0 - DCOS(pi * (DBLE(k1)-1.d0) / (DBLE(ns_part1) * 2.d0))) 
            ds(k1) = xs(k1) - xs(k1-1)
            !write(*,*)k1, xs(k1)
        enddo
     
        !ns_part2 = 2 * NINT(0.5d0*(chainlen - xs_crit)/ds(ns_part1+1))
        !do k1 = ns_part1+2, ns_part1+ns_part2+2
        do k1 = ns_part1+2, ns+1
            xs(k1) = xs(k1-1) + ds(ns_part1+1) 
            ds(k1) = xs(k1) - xs(k1-1)
        enddo

        !do k1 = 1, ns+1
        !   write(123,'(I10,2F10.5)') k1, xs(k1), ds(k1)
        !enddo
    endif
endif

call get_interp_quad_coeff(ds, ns, coeff)

return
!------------------------------------------------------------------------------------------------------!
end subroutine init_chain_contour
