subroutine export_computes(iter, convergence)
!-----------------------------------------------------------------------------------------------------------!
use arrays,       only: phia_mx, phia_gr, wa, wa_new, wa_mix, qmx_final, qmx_interp_mm, qmx_interp_mg, qgr_interp, ds_gr_ed, &
                     &  xs_gr_ed, xs_gr_conv, phia_gr_indiv, coeff_gr_conv, qgr_final, volnp
use hist,         only: planar_cell_of_np, sph_cell_of_np, nbin, lbin, dist_from_face, dist_from_np, cell_vol_planar, cell_vol_sph
use delta,        only: gp_init_value, gpid, num_gpoints
use geometry,     only: numnp, xc, is_dir_face, node_in_q0_face
use parser_vars,  only: gr_exist, mx_exist, n_nanopart_faces, ns_gr_conv, ns_gr_ed, ns_mx_conv, ns_mx_ed,       &
                     &  export_phi_gen_freq, export_field_freq, export_propagators_freq, export_phi_indiv_freq, &
                     &  export_brush_thickness_freq, export_chains_per_area_freq, export_ads_free_freq 
use write_helper, only: adjl, export
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
logical, intent(in) :: convergence

integer, intent(in) :: iter
integer :: mm, nn

character(40) :: file_name
!-----------------------------------------------------------------------------------------------------------!
if (export(export_phi_gen_freq, iter, convergence)) call export_phi_nodal(phia_mx, phia_gr, numnp, xc, volnp)
if (export(export_field_freq, iter, convergence))   call export_field(wa, wa_new, wa_mix)

if (export(export_phi_gen_freq, iter, convergence)) then 
    !chains grafted to planar surfaces
    do mm = 1, 3
        do nn = 1, 2
            if (is_dir_face(mm,nn)) then
                file_name = ""
                write(file_name,'("o.phi_smear_w",I1,"_",I1,".out.txt")') mm, nn
                call export_phi_smear(planar_cell_of_np(:,mm,nn), cell_vol_planar(:,mm,nn), numnp, file_name, phia_mx, phia_gr, volnp, lbin, nbin)
            endif
        enddo
    enddo

    !chains grafted to nanoparticles
    do mm = 1, n_nanopart_faces
        file_name = ""
        write(file_name,'("o.phi_smear_np",I1,".out.txt")') mm
        call export_phi_smear(sph_cell_of_np(mm,:), cell_vol_sph(mm,:), numnp, file_name, phia_mx, phia_gr, volnp, lbin, nbin)
    enddo
endif

if (mx_exist.eq.1) then
    if (export(export_phi_gen_freq, iter, convergence))     call compute_phi_end_middle_nodal(ns_mx_ed, qmx_final, qmx_final, "matrix", numnp, xc)
    if (export(export_propagators_freq, iter, convergence)) call export_q(ns_mx_ed, qmx_final, "mtrx")
endif

if (gr_exist.eq.1) then
    if (export(export_phi_gen_freq, iter, convergence))     call compute_phi_end_middle_nodal(ns_gr_conv, qgr_interp, qmx_interp_mg, "graftd", numnp, xc)
    if (export(export_propagators_freq, iter, convergence)) call export_q(ns_gr_ed, qgr_final, "graf")

    if (export(export_phi_indiv_freq, iter, convergence)) then
        call compute_phi_indiv(numnp, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, wa, num_gpoints, gpid, gp_init_value, phia_gr_indiv)
        call export_phi_indiv(num_gpoints, numnp, xc, phia_gr_indiv)
    endif

    !chains grafted to planar surfaces
    do mm = 1, 3
        do nn = 1, 2
            if (is_dir_face(mm,nn)) then
                if (export(export_brush_thickness_freq, iter, convergence)) then
                    file_name = ""
                    write(file_name,'("o.brush_w",I1,"_",I1,".out.txt")') mm, nn
                    call export_brush(num_gpoints, numnp, phia_gr, phia_gr_indiv, volnp, file_name, dist_from_face(:,mm,nn))
                    file_name = ""
                    write(file_name,'("o.brush99_w",I1,"_",I1,".out.txt")') mm, nn
                    call export_brush99(planar_cell_of_np(:,mm,nn), num_gpoints, numnp, file_name, phia_gr, phia_gr_indiv, volnp, lbin, nbin)
                endif
                if (export(export_chains_per_area_freq, iter, convergence)) then
                    !file_name = ""
                    !write(file_name,'("o.chains_area_w",I1,"_",I1,".out.txt")') mm, nn
                    !call export_chains_area()
                endif
                if (export(export_ads_free_freq, iter, convergence)) then
                endif
            endif
        enddo
    enddo

    !chains grafted to nanoparticles
    do mm = 1, n_nanopart_faces
        if (export(export_brush_thickness_freq, iter, convergence)) then
            file_name = ""
            write(file_name,'("o.brush_np",I1,".out.txt")') mm
            call export_brush(num_gpoints, numnp, phia_gr, phia_gr_indiv, volnp, file_name, dist_from_np(mm,:))
            file_name = ""
            write(file_name,'("o.brush99_np",I1,".out.txt")') mm
            call export_brush99(sph_cell_of_np(mm,:), num_gpoints, numnp, file_name, phia_gr, phia_gr_indiv, volnp, lbin, nbin)
        endif
        if (export(export_chains_per_area_freq, iter, convergence)) then
            !file_name = ""
            !write(file_name,'("o.chains_area_w",I1,"_",I1,".out.txt")') mm, nn
            !call export_chains_area(numnp, node_in_q0_face, sph_cell_of_np(mm,:), ds_gr_ed, num_gpoints, gpid, gp_init_value)
        endif
        if (export(export_ads_free_freq, iter, convergence)) then
        endif
    enddo
#ifdef DEBUG_OUTPUTS
    call export_q(ns_gr_conv, qgr_interp, "grcv")
#endif
endif

!print the restricted partition functions
#ifdef DEBUG_OUTPUTS
call export_q(ns_mx_conv, qmx_interp_mm, "mxcf")
call export_q(ns_gr_conv, qmx_interp_mg, "mxcg")
write(6,'(2X,A40)')adjl("****************************************",40)
#endif
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_computes
