subroutine export_chains_area(numnp, node_in_q0_face, cell_of_np, ds_gr_ed, num_gpoints, gpid, gp_init_value)
!------------------------------------------------------------------------------------------------------!
use parser_vars,  only: ns_gr_ed, mumps_matrix_type
!use hist,         only: nbin
use write_helper, only: adjl
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in)                         :: numnp, num_gpoints
integer, intent(in), dimension(numnp)       :: cell_of_np
integer, intent(in), dimension(num_gpoints) :: gpid
integer                                     :: bin, kk, ii, gnode_id

logical, intent(in), dimension(numnp) :: node_in_q0_face
logical, dimension(numnp)             :: node_in_q0_face_new

real(8), intent(in), dimension(num_gpoints) :: gp_init_value
real(8), intent(in), dimension(ns_gr_ed+1)  :: ds_gr_ed
real(8), dimension(2,numnp)                 :: qgr
real(8), dimension(ns_gr_ed+1,numnp)        :: qgr_final
!------------------------------------------------------------------------------------------------------!
write(6,'(6X,A40)')adjl("*chains per area..",40)
!------------------------------------------------------------------------------------------------------!
do bin = 1, 3
    write(6,*)"bin = ",bin

    node_in_q0_face_new = node_in_q0_face
    do kk = 1, numnp
        if (cell_of_np(kk) .eq. bin) node_in_q0_face_new(kk) = .True.
    enddo

    qgr       = 0.d0
    qgr_final = 0.d0

    do ii = 1, num_gpoints
        gnode_id = gpid(ii)

        qgr(1,gnode_id)       = gp_init_value(ii)
        qgr_final(1,gnode_id) = gp_init_value(ii)
    enddo

    call solver_edwards(ds_gr_ed, ns_gr_ed, mumps_matrix_type, qgr, qgr_final, node_in_q0_face_new)

    !integrate
    !sum_final = 0.d0
     !    sum_Q= 0.d0
!    do mm = 0, nx
!        sum_final = sum_final + coeff_x(mm) * qshape_final(mm,ns) *
!layer_area(mm)
!        sum_Q = sum_Q + coeff_x(mm) * q_final(mm,ns) * layer_area(mm)
!    enddo
!
!    p_cross(kk)  = 1.0 - sum_final / sum_Q

enddo

!sum_phi = 0.0d0
!do kk = 0, nx
!    sum_phi = sum_phi + phi(kk) * coeff_x(kk) * layer_area(kk)
!enddo
!
!do ii=0,nx
!    n_shape(ii) = p_cross(ii) * rho_seg_bulk / chainlen * sum_phi * 1.e-30
!/layer_area(ii)
!enddo
!------------------------------------------------------------------------------------------------------!
end subroutine export_chains_area
