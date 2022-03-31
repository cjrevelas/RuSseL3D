!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_computes(iter, convergence)
!-----------------------------------------------------------------------------------------------------------!
use arrays_mod,       only: phi_mx, phi_gr, phi_gr_indiv, ww, ww_new, ww_mix,                         &
                         &  qmx_final, qgr_final, qmx_interp_mm, qmx_interp_mg, qgr_interp, ds_gr_ed, &
                         &  ds_mx_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, volnp
use hist_mod,         only: nbin, lbin, planar_cell_of_np, sph_cell_of_np, dist_from_face,            &
                         &  dist_from_np, cell_vol_planar, cell_vol_sph
use delta_mod,        only: gp_init_value, gpid, num_gpoints
use geometry_mod,     only: numnp, xc, is_dirichlet_face, node_belongs_to_dirichlet_face
use write_helper_mod, only: adjl, export
use parser_vars_mod,  only: num_of_nanoparticle_faces, mx_exist, gr_exist, ads_distance,              &
                         &  Rg2_per_mon_mx, Rg2_per_mon_gr, chainlen_mx, chainlen_gr,                 &
                         &  ns_mx_ed, ns_mx_conv, ns_gr_ed, ns_gr_conv,                               &
                         &  Rg2_per_mon_mx, Rg2_per_mon_gr, chainlen_mx, chainlen_gr,                 &
                         &  export_phi_gen_freq,                                                      &
                         &  export_field_freq,                                                        &
                         &  export_propagators_freq,                                                  &
                         &  export_phi_indiv_freq,                                                    &
                         &  export_brush_thickness_freq,                                              &
                         &  export_chains_per_area_freq,                                              &
                         &  export_ads_free_freq
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
logical, intent(in)       :: convergence
logical, dimension(numnp) :: adsorbed

integer, intent(in) :: iter
integer             :: kk, mm, nn

character(40) :: file_name
!-----------------------------------------------------------------------------------------------------------!
adsorbed = .false.

if (export(export_phi_gen_freq, iter, convergence)) call export_phi_nodal(phi_mx, phi_gr, numnp, xc, volnp)
if (export(export_field_freq, iter, convergence))   call export_field(ww, ww_new, ww_mix)

if (export(export_phi_gen_freq, iter, convergence)) then
    ! Planar surfaces
    do mm = 1, 3
        do nn = 1, 2
            if (is_dirichlet_face(mm,nn)) then
                file_name = ""
                write(file_name,'("o.phi_smear_w",I1,"_",I1)') mm, nn
                call export_phi_smeared(planar_cell_of_np(:,mm,nn), cell_vol_planar(:,mm,nn), numnp, file_name, phi_mx, phi_gr, volnp, lbin, nbin)
            endif
        enddo
    enddo

    ! Spherical nanoparticles
    do mm = 1, num_of_nanoparticle_faces
        file_name = ""
        write(file_name,'("o.phi_smear_np",I1)') mm
        call export_phi_smeared(sph_cell_of_np(mm,:), cell_vol_sph(mm,:), numnp, file_name, phi_mx, phi_gr, volnp, lbin, nbin)
    enddo
endif

if (mx_exist.eq.1) then
    if (export(export_phi_gen_freq, iter, convergence))     call export_phi_end_middle_nodal(ns_mx_ed, qmx_final, qmx_final, "mx", numnp, xc)
    if (export(export_propagators_freq, iter, convergence)) call export_propagator(ns_mx_ed, qmx_final, "mx")

    ! Planar surfaces
    do mm = 1, 3
        do nn = 1, 2
            if (is_dirichlet_face(mm,nn)) then
                if (export(export_chains_per_area_freq, iter, convergence)) then
                    file_name = ""
                    write(file_name,'("o.chains_area_w",I1,"_",I1)') mm, nn
                    call export_chains_area(node_belongs_to_dirichlet_face, planar_cell_of_np(:,mm,nn), "mx", Rg2_per_mon_mx, chainlen_mx, ns_mx_ed, ds_mx_ed, qmx_final, phi_mx, ww)
                endif
                if (export(export_ads_free_freq, iter, convergence)) then
                    do kk = 1, numnp
                        if (dist_from_face(kk,mm,nn)<ads_distance) adsorbed(kk) = .true.
                    enddo
                endif
            endif
        enddo
    enddo

    ! Spherical nanoparticles
    do mm = 1, num_of_nanoparticle_faces
        if (export(export_chains_per_area_freq, iter, convergence)) then
            file_name = ""
            write(file_name,'("o.chains_area_w",I1,"_",I1)') mm, nn
            call export_chains_area(node_belongs_to_dirichlet_face, sph_cell_of_np(mm,:), "mx", Rg2_per_mon_mx, chainlen_mx, ns_mx_ed, ds_mx_ed, qmx_final, phi_mx, ww)
        endif
        if (export(export_ads_free_freq, iter, convergence)) then
            do kk = 1, numnp
                if (dist_from_np(mm,kk)<ads_distance) adsorbed(kk) = .true.
            enddo
        endif
    enddo

    if (export(export_ads_free_freq, iter, convergence)) call export_ads_free(node_belongs_to_dirichlet_face, adsorbed)

#ifdef DEBUG_OUTPUTS
    call export_propagator(ns_mx_conv, qmx_interp_mm, "mm")
    call export_propagator(ns_gr_conv, qmx_interp_mg, "mg")
    write(6,'(2X,A40)')adjl("****************************************",40)
#endif
endif

if (gr_exist.eq.1) then
    if (export(export_phi_gen_freq, iter, convergence))     call export_phi_end_middle_nodal(ns_gr_conv, qgr_interp, qmx_interp_mg, "gr", numnp, xc)
    if (export(export_propagators_freq, iter, convergence)) call export_propagator(ns_gr_ed, qgr_final, "gr")

    if (export(export_phi_indiv_freq, iter, convergence)) then
        call compute_phi_indiv(numnp, qmx_interp_mg, ds_gr_ed, xs_gr_ed, xs_gr_conv, coeff_gr_conv, ww, num_gpoints, gpid, gp_init_value, phi_gr_indiv)
        call export_phi_indiv(num_gpoints, numnp, xc, phi_gr_indiv)
    endif

    ! Planar surfaces
    do mm = 1, 3
        do nn = 1, 2
            if (is_dirichlet_face(mm,nn)) then
                if (export(export_brush_thickness_freq, iter, convergence)) then
                    file_name = ""
                    write(file_name,'("o.brush_w",I1,"_",I1)') mm, nn
                    call export_brush(num_gpoints, numnp, phi_gr, phi_gr_indiv, volnp, file_name, dist_from_face(:,mm,nn))
                    file_name = ""
                    write(file_name,'("o.brush99_w",I1,"_",I1)') mm, nn
                    call export_brush99(planar_cell_of_np(:,mm,nn), num_gpoints, numnp, file_name, phi_gr, phi_gr_indiv, volnp, lbin, nbin)
                endif
                if (export(export_chains_per_area_freq, iter, convergence)) then
                    file_name = ""
                    write(file_name,'("o.chains_area_w",I1,"_",I1)') mm, nn
                    call export_chains_area(node_belongs_to_dirichlet_face, planar_cell_of_np(:,mm,nn), "gr", Rg2_per_mon_gr, chainlen_gr, ns_gr_ed, ds_gr_ed, qgr_final, phi_gr, ww)
                endif
            endif
        enddo
    enddo

    ! Spherical nanoparticles
    do mm = 1, num_of_nanoparticle_faces
        if (export(export_brush_thickness_freq, iter, convergence)) then
            file_name = ""
            write(file_name,'("o.brush_np",I1)') mm
            call export_brush(num_gpoints, numnp, phi_gr, phi_gr_indiv, volnp, file_name, dist_from_np(mm,:))
            file_name = ""
            write(file_name,'("o.brush99_np",I1)') mm
            call export_brush99(sph_cell_of_np(mm,:), num_gpoints, numnp, file_name, phi_gr, phi_gr_indiv, volnp, lbin, nbin)
        endif
        if (export(export_chains_per_area_freq, iter, convergence)) then
            file_name = ""
            write(file_name,'("o.chains_area_w",I1,"_",I1)') mm, nn
            call export_chains_area(node_belongs_to_dirichlet_face, sph_cell_of_np(mm,:), "gr", Rg2_per_mon_gr, chainlen_gr, ns_gr_ed, ds_gr_ed, qgr_final, phi_gr, ww)
        endif
    enddo
#ifdef DEBUG_OUTPUTS
    call export_propagator(ns_gr_conv, qgr_interp, "gg")
#endif
endif

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine export_computes
