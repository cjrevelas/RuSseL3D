subroutine export_ads_free(node_in_q0_face, adsorbed)
!-----------------------------------------------------------------------------------------------------------------------!
use parser_vars,  only: mumps_matrix_type, Rg2_per_mon_mx, chainlen_mx, ns_mx_ed, ns_mx_conv 
use write_helper, only: adjl
use arrays,       only: ds_mx_ed, xs_mx_ed, xs_mx_conv, coeff_mx_conv, qmx_interp_mm, phia_mx, wa
use geometry,     only: numnp, xc
!-----------------------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------------------!
integer                                :: ii, kk

logical, intent(in), dimension(numnp)  :: node_in_q0_face, adsorbed
logical, dimension(numnp)              :: node_in_q0_face_new

real(8), dimension(2,numnp)            :: qfree
real(8), dimension(ns_mx_ed+1,numnp)   :: qfree_final
real(8), dimension(ns_mx_conv+1,numnp) :: qfree_interp, qads
real(8), dimension(numnp)              :: phi_free, phi_ads, phi_loop, phi_tail
!-----------------------------------------------------------------------------------------------------------------------!
write(6,'(2X,A40)')adjl("Exporting ads vs free density profiles.",40)

!assembly
call fem_matrix_assemble(Rg2_per_mon_mx, wa)

!initial and boundary conditions
qfree               = 0.d0
qfree_final         = 0.d0
qfree(1,:)          = 1.d0
qfree_final(1,:)    = 1.d0
node_in_q0_face_new = node_in_q0_face

do kk = 1, numnp
    if (adsorbed(kk)) then
        qfree(1,kk)             = 0.d0
        qfree_final(1,kk)       = 0.d0
        node_in_q0_face_new(kk) = .true.
    endif
enddo

!solution
call solver_edwards(ds_mx_ed, ns_mx_ed, mumps_matrix_type, qfree, qfree_final, node_in_q0_face_new)

!contour interpolation
do ii = 1, numnp
    call interp_linear(1, ns_mx_ed+1, xs_mx_ed, qfree_final(:,ii), ns_mx_conv+1, xs_mx_conv, qfree_interp(:,ii))
enddo

do ii = 1, ns_mx_conv+1
    do kk = 1, numnp
        qads(ii,kk) = qmx_interp_mm(ii,kk) - qfree_interp(ii,kk)
    enddo
enddo

!determine profiles
call contour_convolution(numnp, chainlen_mx, ns_mx_conv, coeff_mx_conv, qfree_interp, qfree_interp, phi_free)
call contour_convolution(numnp, chainlen_mx, ns_mx_conv, coeff_mx_conv, qads, qads, phi_loop)
call contour_convolution(numnp, chainlen_mx, ns_mx_conv, coeff_mx_conv, qfree_interp, qads, phi_tail)

write(123,'(9A19)') "kk", 'x', 'y', 'z', "phi_free", "phi_ads", "phi_loop", "phi_tail", "phi_mx"
do kk = 1, numnp
    phi_ads(kk) = phia_mx(kk) - phi_free(kk)
    write(123,'(I19,8(E19.9E3))') kk, xc(1,kk), xc(2,kk), xc(3,kk), phi_free(kk), phi_ads(kk), phi_loop(kk), phi_tail(kk), phia_mx(kk)
enddo

return
!-----------------------------------------------------------------------------------------------------------------------!
end subroutine export_ads_free