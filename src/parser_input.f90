!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine parser_input()
!--------------------------------------------------------------------------------!
use parser_vars
use write_helper
use error_handing
use flags
use eos
use iofiles
!--------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------!
character(100) :: line
character(1)   :: axis

integer :: Reason, ii, id, face1, face2

logical :: log_mesh_filename               = .false.
logical :: log_gp_filename                 = .false.
logical :: log_field_in_filename           = .false.
logical :: log_temperature                 = .false.
logical :: log_pressure                    = .false.
logical :: log_mass_density                = .false.
logical :: log_monomer_mass                = .false.
logical :: log_number_of_iterations        = .false.
logical :: log_set_initial_iteration       = .false.
logical :: log_maximum_error               = .false.
logical :: log_maximum_energy_error_tol    = .false.
logical :: log_field_init_scheme           = .false.
logical :: log_fraction_of_new_field       = .false.
logical :: log_mx_exist                    = .false.
logical :: log_Rg2_per_mon_mx              = .false.
logical :: log_ds_ave_mx                   = .false.
logical :: log_chainlen_mx                 = .false.
logical :: log_contour_discr_mx            = .false.
logical :: log_ads_distance                = .false.
logical :: log_gr_exist                    = .false.
logical :: log_Rg2_per_mon_gr              = .false.
logical :: log_chainlen_gr                 = .false.
logical :: log_ds_ave_gr                   = .false.
logical :: log_contour_discr_gr            = .false.
logical :: log_grafted_ic_from_delta       = .false.
logical :: log_num_gr_chains_tol           = .false.
logical :: log_r_gpoint                    = .false.
logical :: log_calc_delta_every            = .false.
logical :: log_n_dirichlet_faces           = .false.
logical :: log_n_nanopart_faces            = .false.
logical :: log_periodicity                 = .false.
logical :: log_sigma_polymer               = .false.
logical :: log_Hamaker_constant_of_polymer = .false.
logical :: log_wall_distance               = .false.
logical :: log_eos_type                    = .false.
logical :: log_eos_coeffs                  = .false.
logical :: log_influence_param             = .false.
logical :: log_profile_dimension           = .false.
logical :: log_mumps_matrix_type           = .false.
logical :: log_bin_thickness               = .false.
logical :: log_export_phi_gen_freq         = .false.
logical :: log_export_phi_indiv_freq       = .false.
logical :: log_export_field_freq           = .false.
logical :: log_export_field_bin_freq       = .false.
logical :: log_export_propagators_freq     = .false.
logical :: log_export_brush_thickness_freq = .false.
logical :: log_export_chains_per_area_freq = .false.
logical :: log_export_ads_free_freq        = .false.
!--------------------------------------------------------------------------------!
! Initialize periodicity arrays
do ii = 1, 3
    periodic_axis_id = .false.
enddo
do ii = 1, 6
    periodic_face_id = -1
enddo

! Parse input file to retrieve simulation parameters
call GET_COMMAND_ARGUMENT(1,input_filename)

if (input_filename == '') input_filename = "./in.input"
inquire(file=input_filename, exist=file_exists)

if (file_exists) then
    open(unit=256, file = input_filename)
else
    write(ERROR_MESSAGE,'("File ",A15," does not exist!")')input_filename
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

do
    read(256,'(A100)',iostat=Reason) line

    if (Reason > 0)  then
        write(*,*)"Something went wrong!"
    elseif (Reason < 0) then
        exit
    else
        if (INDEX(line,"# mesh input file") > 0) then
            read(line,'(A40)') mesh_filename
            log_mesh_filename = .true.
        elseif (INDEX(line,"# gpoints input file") > 0) then
            read(line,'(A40)') gp_filename
            log_gp_filename = .true.
        elseif (INDEX(line,"# field input file") > 0) then
            read(line,*) field_in_filename
            log_field_in_filename = .true.
        elseif (INDEX(line,"# temp") > 0) then
            read(line,*) temp
            log_temperature = .true.
        elseif (INDEX(line," pres") > 0) then
            read(line,*) pres
            log_pressure = .true.
        elseif (INDEX(line,"# mass den") > 0) then
            read(line,*) massden
            log_mass_density = .true.
        elseif (INDEX(line,"# mon mass") > 0) then
            read(line,*) mon_mass
            log_monomer_mass = .true.
        elseif (INDEX(line,"# Rg2/mon matrix") > 0) then
            read(line,*) Rg2_per_mon_mx
            log_Rg2_per_mon_mx = .true.
        elseif (INDEX(line,"# Rg2/mon grafted") > 0) then
            read(line,*) Rg2_per_mon_gr
            log_Rg2_per_mon_gr = .true.
        elseif (INDEX(line,"# chain length matrix") > 0) then
            read(line,*) chainlen_mx
            log_chainlen_mx = .true.
        elseif (INDEX(line,"# chain length grafted") > 0) then
            read(line,*) chainlen_gr
            log_chainlen_gr = .true.
        elseif (INDEX(line,"# num gr chains tol") > 0) then
            read(line,*) num_gr_chains_tol
            log_num_gr_chains_tol = .true.
        elseif (INDEX(line,"# pol sigma") > 0) then
            read(line,*) sigma_pol
            log_sigma_polymer = .true.
        elseif (INDEX(line,"# pol ham") > 0) then
            read(line,*) A_pol
            log_Hamaker_constant_of_polymer = .true.
        elseif (INDEX(line,"# wall dist") > 0) then
            read(line,*) wall_distance
            log_wall_distance = .true.
        elseif (INDEX(line,"# fraction") > 0) then
            read(line,*) frac
            log_fraction_of_new_field = .true.
        elseif (INDEX(line,"# max error") > 0) then
            read(line,*) max_error_tol
            log_maximum_error = .true.
        elseif (INDEX(line,"# max energy error") > 0) then
            read(line,*) free_energy_error_tol
            log_maximum_energy_error_tol = .true.
        elseif (INDEX(line,"# init iter") > 0) then
            read(line,*) init_iter
            log_set_initial_iteration= .true.
        elseif (INDEX(line,"# num iter") > 0) then
            read(line,*) iterations
            log_number_of_iterations = .true.
        elseif (INDEX(line,"# init field") > 0) then
            read(line,*) field_init_scheme
            log_field_init_scheme = .true.
        elseif (INDEX(line,"# bin thickness") > 0) then
            read(line,*) bin_thickness
            log_bin_thickness = .true.
        elseif (INDEX(line,"# export dens profs") > 0) then
            read(line,*) export_phi_gen_freq
            log_export_phi_gen_freq = .true.
        elseif (INDEX(line,"# export indiv dens profs") > 0) then
            read(line,*) export_phi_indiv_freq
            log_export_phi_indiv_freq = .true.
        elseif (INDEX(line,"# export field") > 0) then
            read(line,*) export_field_freq
            log_export_field_freq = .true.
        elseif (INDEX(line,"# export binary field") > 0) then
            read(line,*) export_field_bin_freq
            log_export_field_bin_freq = .true.
        elseif (INDEX(line,"# export propagators") > 0) then
            read(line,*) export_propagators_freq
            log_export_propagators_freq = .true.
        elseif (INDEX(line,"# export brush thickness") > 0) then
            read(line,*) export_brush_thickness_freq
            log_export_brush_thickness_freq = .true.
        elseif (INDEX(line,"# export chains per area profs") > 0) then
            read(line,*) export_chains_per_area_freq
            log_export_chains_per_area_freq = .true.
        elseif (INDEX(line,"# export ads vs free profs") > 0) then
            read(line,*) export_ads_free_freq
            log_export_ads_free_freq = .true.
        elseif (INDEX(line,"# prof dim") > 0) then
            read(line,*) prof_dim
            log_profile_dimension = .true.
        elseif (INDEX(line,"# use matrix") > 0) then
            read(line,*) mx_exist
            log_mx_exist = .true.
        elseif (INDEX(line,"# use grafted") > 0) then
            read(line,*) gr_exist
            log_gr_exist = .true.
        elseif (INDEX(line,"# contour step matrix") > 0) then
            read(line,*) ds_ave_mx_ed, ds_ave_mx_conv
            log_ds_ave_mx = .true.
        elseif (INDEX(line,"# discr scheme matrix") > 0) then
            read(line,*) contour_discr_mx
            if (contour_discr_mx.eq.contour_hybrid) read(256,*) xs_crit_mx
            log_contour_discr_mx = .true.
        elseif (INDEX(line,"# ads distance") > 0) then
            read(line,*) ads_distance
            log_ads_distance = .true.
        elseif (INDEX(line,"# contour step grafted") > 0) then
            read(line,*) ds_ave_gr_ed, ds_ave_gr_conv
            log_ds_ave_gr = .true.
        elseif (INDEX(line,"# gp dist from solid") > 0) then
            read(line,*) r_gpoint
            log_r_gpoint = .true.
        elseif (INDEX(line,"# discr scheme grafted") > 0) then
            read(line,*) contour_discr_gr
            if (contour_discr_gr.eq.contour_hybrid) read(256,*) xs_crit_gr
            log_contour_discr_gr = .true.
        elseif (INDEX(line,"# mumps matrix") > 0) then
            read(line,*) mumps_matrix_type
            log_mumps_matrix_type = .true.
        elseif (INDEX(line,"# eos type") > 0) then
            read(line,*) eos_type
            log_eos_type = .true.
        elseif (INDEX(line,"# eos coeffs") > 0) then
            if (eos_type.eq.eos_helfand) then
                read(line,*) hlf_kappa_T
            elseif (eos_type.eq.eos_sl)  then
                read(line,*) rho_star, T_star, P_star
            endif
            log_eos_coeffs = .true.
        elseif (INDEX(line, "# eos infl param") > 0) then
            read(line,*) k_gr_tilde
            log_influence_param = .true.
        elseif (INDEX(line,"# calc delta") > 0) then
            read(line,*) grafted_ic_from_delta
            log_grafted_ic_from_delta = .true.
        elseif (INDEX(line,"# freq delta") > 0) then
            read(line,*) calc_delta_every
            log_calc_delta_every = .true.
        elseif (INDEX(line,"# num faces") > 0) then
            read(line,*) n_dirichlet_faces
            if (n_dirichlet_faces > 0) then
                allocate(ids_dirichlet_faces(n_dirichlet_faces))
                allocate(sigma_plate(n_dirichlet_faces))
                allocate(A_plate(n_dirichlet_faces))
                do ii = 1, n_dirichlet_faces
                    read(256,*) id, sigma_plate(ii), A_plate(ii)
                    ids_dirichlet_faces(ii) = id
                enddo
                A_plate = A_plate * 1e-20
                log_n_dirichlet_faces = .true.
            endif
        elseif (INDEX(line,"# num nanop") > 0) then
            read(line,*) n_nanopart_faces
            if (n_nanopart_faces > 0) then
                allocate(ids_nanopart_faces(n_nanopart_faces))
                allocate(center_np(3,n_nanopart_faces))
                allocate(radius_np_eff(n_nanopart_faces))
                allocate(sigma_np(n_nanopart_faces))
                allocate(A_np(n_nanopart_faces))
                do ii = 1, n_nanopart_faces
                    read(256,*) id, radius_np_eff(ii), center_np(1,ii), center_np(2,ii), center_np(3,ii), &
                                                          & sigma_np(ii), A_np(ii)
                    ids_nanopart_faces(ii) = id
                enddo
                A_np = A_np * 1e-20
                log_n_nanopart_faces = .true.
            endif
        elseif (INDEX(line,"# periodicity") > 0) then
            read(line,*) periodicity
            do ii = 1, periodicity
                read(256,*) axis, face1, face2
                if (axis=='x') then
                    periodic_axis_id(1) = .true.
                    periodic_face_id(1) = face1
                    periodic_face_id(2) = face2
                elseif (axis=='y') then
                    periodic_axis_id(2) = .true.
                    periodic_face_id(3) = face1
                    periodic_face_id(4) = face2
                elseif (axis=='z') then
                    periodic_axis_id(3) = .true.
                    periodic_face_id(5) = face1
                    periodic_face_id(6) = face2
                else
                    write(ERROR_MESSAGE,'("Periodicity axis is not valid. Must be x, y or z.")')
                    call exit_with_error(1,1,1,ERROR_MESSAGE)
                endif
            enddo
            log_periodicity = .true.
        endif
   endif
enddo

close(256)

! Check input parameters
write(iow,'(A85)')adjl('-----------------------------------SYSTEM PARAMETERS---------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------SYSTEM PARAMETERS---------------------------------',85)


if (log_temperature) then
    if (Temp>0) then
        write(iow,'(3X,A40,E16.9,A4)')adjl("Temperature:",40),Temp," [K]"
        write(6  ,'(3X,A40,E16.9,A4)')adjl("Temperature:",40),Temp," [K]"
        beta = 1.d0 / (boltz_const_Joule_K * Temp)
    else
        write(ERROR_MESSAGE,'("Temperature is negative: ",E16.9," K")') temp
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Temperature was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif


if (log_pressure) then
    if (pres>=0) then
        write(iow,'(3X,A40,E16.9,A6)')adjl("Pressure:",40),pres," [atm]"
        write(*  ,'(3X,A40,E16.9,A6)')adjl("Pressure:",40),pres," [atm]"
        pres = pres * atm_to_pa
    else
        write(ERROR_MESSAGE,'("Pressure is negative: ",E16.9, " atm")') pres
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    pres = 0.d0
    write(iow,'(A40)') 'Pressure was set to 0 atm'
    write(*  ,'(A40)') 'Pressure was set to 0 atm'
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('----------------------------------POLYMER PROPERTIES---------------------------------',85)
write(*  ,'(A85)')adjl('----------------------------------POLYMER PROPERTIES---------------------------------',85)


if (log_mass_density) then
    if (massden>0) then
        write(iow,'(3X,A40,E16.9,A8)')adjl("Mass density:",40),massden," [g/cm3]"
        write(6  ,'(3X,A40,E16.9,A8)')adjl("Mass density:",40),massden," [g/cm3]"
    else
        write(ERROR_MESSAGE,'("Mass density is negative: ",E16.9," g/cm3")') massden
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Mass density was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif


if (log_monomer_mass) then
    if (mon_mass>0) then
        write(iow,'(3X,A40,E16.9,A8)')adjl("Monomer mass:",40),mon_mass,"[g/mol]"
        write(6  ,'(3X,A40,E16.9,A8)')adjl("Monomer mass:",40),mon_mass,"[g/mol]"
    else
        write(ERROR_MESSAGE,'("Monomer mass is negative: ",E16.9)') mon_mass
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Monomer mass not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif


if (mx_exist.eq.1) then
    if (log_chainlen_mx) then
        if (chainlen_mx>0) then

            chainlen_mx_max = chainlen_mx
        else
            write(ERROR_MESSAGE,'("Chain length of matrix chains is negative: ",E16.9)') chainlen_mx
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Chain length of matrix chains was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
endif


if (log_gr_exist) then
    if (gr_exist.eq.1) then
        write(iow,*)
        write(*,*)
        write(iow,'(A85)')adjl('-----------------------------------GRAFTED CHAINS-----------------------------------',85)
        write(*  ,'(A85)')adjl('-----------------------------------GRAFTED CHAINS-----------------------------------',85)
    else
        gr_exist = 0
    endif
else
    continue
endif


if (gr_exist.eq.1) then
    if (log_gp_filename) then
        inquire(file=gp_filename, exist=file_exists)
        if (.not.file_exists) then
            write(ERROR_MESSAGE,'("Grafting points file ",A16," does not exist!")')gp_filename
            call exit_with_error(1,1,1,ERROR_MESSAGE)
            STOP
        endif
        write(iow,'(3X,A40,A16)')adjl("Reading grafting points from file:",40),gp_filename
        write(6  ,'(3X,A40,A16)')adjl("Reading grafting points from file:",40),gp_filename
    else
        gp_filename = "./in.gnodes"
        write(iow,'(3X,A40)')adjl("Grafting points input file not specified.",40)
        write(iow,'(3X,A40,A16)')adjl("Reading default grafting points file:",40),gp_filename
        write(6  ,'(3X,A40)')adjl("Grafting points input file not specified.",40)
        write(6  ,'(3X,A40,A16)')adjl("Reading default grafting points file:",40),gp_filename

        inquire(file=gp_filename, exist=file_exists)
        if (.not.file_exists) then
            write(ERROR_MESSAGE,'("Default grafting points file ",A16," does not exist!")')gp_filename
            call exit_with_error(1,1,1,ERROR_MESSAGE)
            STOP
        endif
    endif

    if (log_Rg2_per_mon_gr) then
        if (Rg2_per_mon_gr>0) then
            write(iow,'(3X,A40,E16.9,A13)')adjl("Rg2 per grafted monomer:",40),Rg2_per_mon_gr,"[Angstrom^2]"
            write(6  ,'(3X,A40,E16.9,A13)')adjl("Rg2 per grafted monomer:",40),Rg2_per_mon_gr,"[Angstrom^2]"
        else
            write(ERROR_MESSAGE,'("Rg2 per grafted monomer is negative: ",E16.9)') Rg2_per_mon_gr
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Rg2 per matrix monomer was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_chainlen_gr) then
        if (chainlen_gr>0) then
            write(iow,'(3X,A40,E16.9,A11)')adjl("Chain length of grafted chains:",40),chainlen_gr,"[monomers]"
            write(iow,'(3X,A40,E16.9,A11)')adjl("Radius of gyration of grafted chains:",40),&
                                                              & DSQRT(Rg2_per_mon_gr*chainlen_gr),"[Angstrom]"
            write(6  ,'(3X,A40,E16.9,A11)')adjl("Chain length of grafted chains:",40),chainlen_gr,"[monomers]"
            write(6  ,'(3X,A40,E16.9,A11)')adjl("Radius of gyration of grafted chains:",40),&
                                                              & DSQRT(Rg2_per_mon_gr*chainlen_gr),"[Angstrom]"
            chainlen_mx_max = MAX(chainlen_mx, chainlen_gr)
        else
            write(ERROR_MESSAGE,'("Chain length of grafted chains is negative: ",E16.9)') chainlen_gr
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Chain length of grafted chains was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_ds_ave_gr) then
        if (ds_ave_gr_ed>0 .and. ds_ave_gr_conv>0) then
            ns_gr_ed   = 2 * NINT(0.5d0 * chainlen_gr / ds_ave_gr_ed)
            ns_gr_conv = 2 * NINT(0.5d0 * chainlen_gr / ds_ave_gr_conv)
            write(iow,'(3X,A40,I9,I7)')adjl("Number of grafted segments:",40), ns_gr_ed, ns_gr_conv
            write(6  ,'(3X,A40,I9,I7)')adjl("Number of grafted segments:",40), ns_gr_ed, ns_gr_conv

            if (MOD(ns_gr_ed,2).ne.0 .or. MOD(ns_gr_conv,2).ne.0) then
                write(ERROR_MESSAGE,'("ns_grafted is not an even number: ",I16,I16)') ns_gr_ed, ns_gr_conv
                call exit_with_error(1,1,1,ERROR_MESSAGE)
            endif
        else
            write(ERROR_MESSAGE,'("Contour step of grafted chains is negative: ",E16.9,E16.9)') ds_ave_gr_ed, ds_ave_gr_conv
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Contour step of grafted chains was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_contour_discr_gr) then
        if (contour_discr_gr.eq.contour_uniform) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "uniform"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "uniform"
        elseif (contour_discr_gr.eq.contour_symm) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "symm"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "symm"
        elseif (contour_discr_gr.eq.contour_hybrid) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "hybrid"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "hybrid"
            if (xs_crit_gr < 0) then
                write(ERROR_MESSAGE,'("Critical contour point of grafted chains is negative: ",E16.9)') xs_crit_gr
                call exit_with_error(1,1,1,ERROR_MESSAGE)
            elseif (xs_crit_gr>0) then
                write(iow,'(3X,A40,F16.9)')adjl("Crit contour point of grafted chains:",40), xs_crit_gr
                write(6  ,'(3X,A40,F16.9)')adjl("Crit contour point of grafted chains:",40), xs_crit_gr
            else
                ERROR_MESSAGE="Critical contour point of grafted chains was not detected."
                call exit_with_error(1,1,1,ERROR_MESSAGE)
            endif

            call get_contour_step(ds_ave_gr_ed, xs_crit_gr, chainlen_gr, ns_gr_ed)

        elseif (contour_discr_gr.eq.contour_asymm) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "asymm"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "asymm"
        else
            write(ERROR_MESSAGE,'("Not valid Edwards contour scheme of grafted chains: ",I5)') contour_discr_gr
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Edwards contour scheme of grafted chains was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
endif


if (log_mx_exist) then
    if (mx_exist.eq.1) then
        write(iow,*)
        write(*,*)
        write(iow,'(A85)')adjl('-------------------------------------MATRIX CHAINS-----------------------------------',85)
        write(*  ,'(A85)')adjl('-------------------------------------MATRIX CHAINS-----------------------------------',85)
    else
        mx_exist = 0
    endif
else
    continue
endif


if (mx_exist.eq.1) then
    if (log_Rg2_per_mon_mx) then
        if (Rg2_per_mon_mx>0) then
            write(iow,'(3X,A40,E16.9,A13)')adjl("Rg2 per matrix monomer:",40),Rg2_per_mon_mx,"[Angstrom^2]"
            write(6  ,'(3X,A40,E16.9,A13)')adjl("Rg2 per matrix monomer:",40),Rg2_per_mon_mx,"[Angstrom^2]"
        else
            write(ERROR_MESSAGE,'("Rg2 per matrix monomer is negative: ",E16.9)') Rg2_per_mon_mx
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Rg2 per matrix monomer was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    write(iow,'(3X,A40,E16.9,A11)')adjl("Chain length of matrix chains:",40),chainlen_mx,"[monomers]"
    write(iow,'(3X,A40,E16.9,A11)')adjl("Radius of gyration of matrix chains:",40),&
                                                             & DSQRT(Rg2_per_mon_mx*chainlen_mx),"[Angstrom]"
    write(6  ,'(3X,A40,E16.9,A11)')adjl("Chain length of matrix chains:",40),chainlen_mx,"[monomers]"
    write(6  ,'(3X,A40,E16.9,A11)')adjl("Radius of gyration of matrix chains:",40),&
                                                             & DSQRT(Rg2_per_mon_mx*chainlen_mx),"[Angstrom]"
    if (log_ds_ave_mx) then
        if (ds_ave_mx_ed>0 .and. ds_ave_mx_conv>0) then
            ns_mx_ed   = 2 * NINT(0.5d0 * chainlen_mx_max / ds_ave_mx_ed)
            ns_mx_conv = 2 * NINT(0.5d0 * chainlen_mx     / ds_ave_mx_conv)

            write(iow,'(3X,A40,I9,I7)')adjl("Number of matrix segments:",40), ns_mx_ed, ns_mx_conv
            write(6  ,'(3X,A40,I9,I7)')adjl("Number of matrix segments:",40), ns_mx_ed, ns_mx_conv

            if (MOD(ns_mx_ed,2).ne.0 .or. MOD(ns_mx_conv,2).ne.0) then
                write(ERROR_MESSAGE,'("ns_matrix is not an even number: ",I16,I16)') ns_mx_ed, ns_mx_conv
                call exit_with_error(1,1,1,ERROR_MESSAGE)
            endif
        else
            write(ERROR_MESSAGE,'("Contour step of matrix chains is negative: ",E16.9,E16.9)') ds_ave_mx_ed, ds_ave_mx_conv
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Contour step of matrix chains was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_contour_discr_mx) then
        if (contour_discr_mx.eq.contour_uniform) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "uniform"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "uniform"
        elseif (contour_discr_mx.eq.contour_symm) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "symm"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "symm"
        elseif (contour_discr_mx.eq.contour_hybrid) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "hybrid"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "hybrid"
            if (xs_crit_mx < 0) then
                write(ERROR_MESSAGE,'("Critical contour point of matrix chains is negative: ",E16.9)') xs_crit_mx
                call exit_with_error(1,1,1,ERROR_MESSAGE)
            elseif (xs_crit_mx > 0) then
                write(iow,'(3X,A40,F16.9)')adjl("Critical contour point of matrix chains:",40), xs_crit_mx
                write(6  ,'(3X,A40,F16.9)')adjl("Critical contour point of matrix chains:",40), xs_crit_mx
            else
                ERROR_MESSAGE="Critical contour point of matrix chains was not detected."
                call exit_with_error(1,1,1,ERROR_MESSAGE)
            endif

            call get_contour_step(ds_ave_mx_ed, xs_crit_mx, chainlen_mx_max, ns_mx_ed)

        elseif (contour_discr_mx.eq.contour_asymm) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "asymm"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "asymm"
        else
            write(ERROR_MESSAGE,'("Not valid Edwards contour scheme of matrix chains: ",I5)') contour_discr_mx
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Edwards contour scheme of matrix chains was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_ads_distance) then
        if (ads_distance.ge.0.d0) then
            write(iow,'(3X,A40,E16.9,A11)')adjl("Adsorption distance for chain segments:",40),ads_distance,"[Angstrom]"
            write(6  ,'(3X,A40,E16.9,A11)')adjl("Adsorption distance for chain segments:",40),ads_distance,"[Angstrom]"
        else
            write(ERROR_MESSAGE,'("Adsorption distance for chain segments is negative:",E16.9)')ads_distance
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ads_distance = 1.2d1
        write(iow,'(3X,A40)')adjl("Adsorption distance not found.",40)
        write(iow,'(3X,A40,E16.9,A11)')adjl("It was set to the default value:",40),ads_distance,"[Angstrom]"
        write(6  ,'(3X,A40)')adjl("Adsorption distance not found.",40)
        write(6  ,'(3X,A40,E16.9,A11)')adjl("It was set to the default value:",40),ads_distance,"[Angstrom]"
    endif
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('--------------------------------SIMULATION PARAMETERS---------------------------------',85)
write(*  ,'(A85)')adjl('--------------------------------SIMULATION PARAMETERS---------------------------------',85)


if (log_mx_exist.and.mx_exist.ge.1) then
    mx_exist = 1
    write(iow,'(3X,A40,I9)')adjl("The system includes matrix chains:",40),mx_exist
    write(6  ,'(3X,A40,I9)')adjl("The system includes matrix chains:",40),mx_exist
else
    mx_exist = 0
    write(iow,'(3X,A40,I9)')adjl("System does not include matrix chains:",40),mx_exist
    write(6  ,'(3X,A40,I9)')adjl("System does not include matrrix chains:",40),mx_exist
endif


if (log_gr_exist.and.gr_exist.ge.1) then
    gr_exist = 1
    write(iow,'(3X,A40,I9)')adjl("The system includes grafted chains:",40),gr_exist
    write(6  ,'(3X,A40,I9)')adjl("The system includes grafted chains:",40),gr_exist

    if (log_grafted_ic_from_delta) then
        write(iow,'(3X,A40,I9)')adjl("Grafted ic from delta:",40),grafted_ic_from_delta
        write(6  ,'(3X,A40,I9)')adjl("Grafted ic from delta:",40),grafted_ic_from_delta
        if (grafted_ic_from_delta.ge.1) then
            if (log_calc_delta_every) then
                if (calc_delta_every.gt.0) then
                    write(iow,'(3X,A40,I9)')adjl("Delta is calculated every:",40),calc_delta_every
                    write(6  ,'(3X,A40,I9)')adjl("Delta is calculated every:",40),calc_delta_every
                else
                    write(iow,'(3X,A40)')adjl("Delta is read from file",40)
                    write(6  ,'(3X,A40)')adjl("Delta is read from file",40)
                endif
            endif
            if (log_num_gr_chains_tol) then
                if (num_gr_chains_tol.ge.0.d0) then
                    write(iow,'(3X,A40,E16.9)')adjl("Number of grafted chains tolerance:",40),num_gr_chains_tol
                    write(6  ,'(3X,A40,E16.9)')adjl("Number of grafted chains tolerance:",40),num_gr_chains_tol
                else
                    write(ERROR_MESSAGE,'("Number of grafted chains tolerance is negative:",E16.9)') num_gr_chains_tol
                    call exit_with_error(1,1,1,ERROR_MESSAGE)
                endif
            else
                num_gr_chains_tol = 0.005d0
                write(iow,'(3X,A40)')adjl("Num grafted chains tolerance not found",40)
                write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),num_gr_chains_tol
                write(6  ,'(3X,A40)')adjl("Num grafted chains tolerance not found",40)
                write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),num_gr_chains_tol
            endif
        endif
        if (grafted_ic_from_delta.lt.1) then
            write(iow,'(3X,A40)')adjl("Initial conditions are read from file.",40)
            write(6  ,'(3X,A40)')adjl("Initial conditions are read from file.",40)
        endif
    else
        write(iow,'(3X,A40)')adjl("Initial conditions are read from file.",40)
        write(6  ,'(3X,A40)')adjl("Initial conditions are read from file.",40)
    endif

    if (log_r_gpoint) then
        if (r_gpoint.gt.0.d0) then
            write(iow,'(3X,A40,E16.9)')adjl("Distance of grafting points from solid:",40),r_gpoint
            write(6  ,'(3X,A40,E16.9)')adjl("Distance of grafting points from solid:",40),r_gpoint
        else
            write(ERROR_MESSAGE,'("Distance of grafting points from solid must have a positive value:",E16.9)') r_gpoint
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        r_gpoint = 0.4d0
        write(iow,'(3X,A40)')adjl("Distance of grafting points not found.",40)
        write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),r_gpoint
        write(6  ,'(3X,A40)')adjl("Distance of grafting points not found.",40)
        write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),r_gpoint
    endif
else
    gr_exist = 0
    write(iow,'(3X,A40,1x,I15)')adjl("System does not include grafted chains:",40),gr_exist
    write(6  ,'(3X,A40,1x,I15)')adjl("System does not include grafted chains:",40),gr_exist
endif


if (log_set_initial_iteration) then
    if (init_iter.lt.0) then
        write(ERROR_MESSAGE,'("Wrong value of initial iteration.",I16)') init_iter
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    init_iter = 0
endif


if (init_iter.eq.0) then
    write(iow,'(3X,A40,I9)')adjl("Fresh simulation starting from iter:",40), init_iter
    write(6  ,'(3X,A40,I9)')adjl("Fresh simulation starting from iter:",40), init_iter
elseif (init_iter.gt.0) then
    write(iow,'(3X,A40,I9)')adjl("Simulation restarting from iter:",40), init_iter
    write(6  ,'(3X,A40,I9)')adjl("Simulation restarting from iter:",40), init_iter
endif


if (log_number_of_iterations) then
    if (iterations.ge.0) then
        write(iow,'(3X,A40,I9)')adjl("Maximum number of iterations:",40),iterations
        write(6  ,'(3X,A40,I9)')adjl("Maximum number of iterations:",40),iterations
    else
        write(ERROR_MESSAGE,'("Maximum number of iterations is negative:",I10)') iterations
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    iterations = 1
    write(iow,'(3X,A40)')adjl("Max number of iter not found.",40)
    write(iow,'(3X,A40,I9)')adjl("It was set to the default value:",40),iterations
    write(6  ,'(3X,A40)')adjl("Max number of iter not found.",40)
    write(6  ,'(3X,A40,I9)')adjl("It was set to the default value:",40),iterations
endif


if (log_bin_thickness) then
    if (bin_thickness.ge.0.d0) then
        write(iow,'(3X,A40,E16.9,A11)')adjl("Bin thickness:",40),bin_thickness,"[Angstrom]"
        write(6  ,'(3X,A40,E16.9,A11)')adjl("Bin thickness:",40),bin_thickness,"[Angstrom]"
    else
        write(ERROR_MESSAGE,'("Bin thickness is negative:",E16.9)')bin_thickness
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    bin_thickness = 0.5d0
    write(iow,'(3X,A40)')adjl("Bin thickness not found.",40)
    write(iow,'(3X,A40,E16.9,A11)')adjl("It was set to the default value:",40),bin_thickness,"[Angstrom]"
    write(6  ,'(3X,A40)')adjl("Bin thickness not found.",40)
    write(6  ,'(3X,A40,E16.9,A11)')adjl("It was set to the default value:",40),bin_thickness,"[Angstrom]"
endif


if (log_profile_dimension) then
    if ( (prof_dim.ge.1).and.(prof_dim.le.3) ) then
        write(iow,'(3X,A40,I9)')adjl("Profile dimension:",40),prof_dim
        write(6  ,'(3X,A40,I9)')adjl("Profile dimension:",40),prof_dim
    else
        write(ERROR_MESSAGE,'("Profile dimension is not between 1 and 3:",I16)') prof_dim
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    prof_dim = 3
    write(iow,'(3X,A40)')adjl("Profile dimension not found.",40)
    write(iow,'(3X,A40,I8)')adjl("It was set to the default value:",40),prof_dim
    write(6  ,'(3X,A40)')adjl("Profile dimension not found.",40)
    write(6  ,'(3X,A40,I8)')adjl("It was set to the default value:",40),prof_dim
endif


if (log_field_init_scheme) then
    if (field_init_scheme==0) then
        write(iow,'(3X,A34)')adjl("Field will be initialized to zero.",40)
        write(6  ,'(3X,A34)')adjl("Field will be initialized to zero.",40)
    elseif (field_init_scheme==1) then
        if (log_gp_filename) then
            inquire(file=field_in_filename, exist=file_exists)
            if (.not.file_exists) then
                write(ERROR_MESSAGE,'("Field input file ",A16," does not exist!")')field_in_filename
                call exit_with_error(1,1,1,ERROR_MESSAGE)
                STOP
            endif
            write(iow,'(A43,A16)')adjl("Field will be read from file:",40),field_in_filename
            write(6  ,'(A43,A16)')adjl("Field will be read from file:",40),field_in_filename
        else
            field_in_filename = "./in.field.bin"
            write(iow,'(3X,A40)')adjl("Field input file not specified.",40)
            write(iow,'(3X,A40,A16)')adjl("Reading default field input file:",40),field_in_filename
            write(6  ,'(3X,A40)')adjl("Field input file not specified.",40)
            write(6  ,'(3X,A40,A16)')adjl("Reading default field input file:",40),field_in_filename

            inquire(file=field_in_filename, exist=file_exists)
            if (.not.file_exists) then
                write(ERROR_MESSAGE,'("Default field input file ",A16," does not exist!")')field_in_filename
                call exit_with_error(1,1,1,ERROR_MESSAGE)
                STOP
            endif
        endif
    elseif (field_init_scheme==2) then
        write(iow,'(3X,A40)')adjl("Field: -kapa at Dir. and 0 elsewhere.",40)
        write(6  ,'(3X,A40)')adjl("Field: -kapa at Dir. and 0 elsewhere.",40)
    else
        write(ERROR_MESSAGE,'("Incorrect field initialization value. Choose between 1-3.",I16)') init_iter
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
endif
if (.not.log_field_init_scheme) then
    write(iow,'(A40)')adjl("Field will be initialized to zero.",40)
    write(6  ,'(A40)')adjl("Field will be initialized to zero.",40)
endif


if (log_mumps_matrix_type) then
    if (mumps_matrix_type == 0) then
        write(iow,'(3X,A40,1X,A14,I1,A1)')adjl("MUMPS matrix type:",40),"nonsymmetric (", mumps_matrix_type,")"
        write(6  ,'(3X,A40,1X,A14,I1,A1)')adjl("MUMPS matrix type:",40),"nonsymmetric (", mumps_matrix_type,")"
    elseif (mumps_matrix_type == 1) then
        write(iow,'(3X,A40,1X,A21,I1,A1)')adjl("MUMPS matrix type:",40),"symmetric pos. def. (", mumps_matrix_type,")"
        write(6  ,'(3X,A40,1X,A21,I1,A1)')adjl("MUMPS matrix type:",40),"symmetric pos. def. (", mumps_matrix_type,")"
    elseif (mumps_matrix_type == 2) then
        write(iow,'(3X,A40,1X,A18,I1,A1)')adjl("MUMPS matrix type:",40),"general symmetric(", mumps_matrix_type,")"
        write(6  ,'(3X,A40,1X,A18,I1,A1)')adjl("MUMPS matrix type:",40),"general symmetric(", mumps_matrix_type,")"
    else
        write(ERROR_MESSAGE,'("Incorrect MUMPS matrix type.",I16)') mumps_matrix_type
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    mumps_matrix_type = 0
    write(iow,'(3X,A40)')adjl("MUMPS matrix type not detected.",40)
    write(iow,'(3X,A40,I10)')adjl("It was set to the nonsymmetric:",40),mumps_matrix_type
    write(6  ,'(3X,A40)')adjl("MUMPS matrix type not detected.",40)
    write(6  ,'(3X,A40,I10)')adjl("It was set to the nonsymmetric:",40),mumps_matrix_type
endif


if (log_fraction_of_new_field) then
    if (frac.ge.0 .and. frac.le.1) then
        write(iow,'(3X,A40,E16.9)')adjl("Initial fraction of new field:",40),frac
        write(6  ,'(3X,A40,E16.9)')adjl("Initial fraction of new field:",40),frac
    else
        write(ERROR_MESSAGE,'("Initial fraction of new field is negative or larger than unity:",E16.9)') frac
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    frac = 1.d0
    write(iow,'(3X,A40)')adjl("No initial fraction of new field.",40)
    write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),frac
    write(6  ,'(3X,A40)')adjl("No initial fraction of new field.",40)
    write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),frac
endif


if (log_maximum_error) then
    if (max_error_tol.ge.0.d0) then
        write(iow,'(3X,A40,E16.9)')adjl("Maximum tolerance error on field:",40),max_error_tol
        write(6  ,'(3X,A40,E16.9)')adjl("Maximum tolerance error on field:",40),max_error_tol
    else
        write(ERROR_MESSAGE,'("Maximum tolerance error on field is negative:",E16.9)') max_error_tol
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    max_error_tol = 0.d0
    write(iow,'(3X,A40)')adjl("Max field error not found",40)
    write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),max_error_tol
    write(6  ,'(3X,A40)')adjl("Max field error not found.",40)
    write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),max_error_tol
endif

if (log_maximum_energy_error_tol) then
    if (free_energy_error_tol.ge.0.d0) then
        write(iow,'(3X,A40,E16.9)')adjl("Maximum tolerance error on energy:",40),free_energy_error_tol
        write(6  ,'(3X,A40,E16.9)')adjl("Maximum tolerance error on energy:",40),free_energy_error_tol
    else
        write(ERROR_MESSAGE,'("Maximum tolerance error is negative:",E16.9)')free_energy_error_tol
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    free_energy_error_tol = 1.d-6
    write(iow,'(3X,A40)')adjl("Max energy error tol not found",40)
    write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),free_energy_error_tol
    write(6  ,'(3X,A40)')adjl("Max energy error tol not found.",40)
    write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),free_energy_error_tol
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-----------------------------------OUTPUT FREQUENCY-----------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------OUTPUT FREQUENCY-----------------------------------',85)


if (log_export_phi_gen_freq) then
    if (export_phi_gen_freq.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export density profiles:",40),export_phi_gen_freq
        write(6  ,'(3X,A40,I9)')adjl("Export density profiles:",40),export_phi_gen_freq
    else
        export_phi_gen_freq = 0
        write(iow,'(3X,A40,I9)')adjl("Export density profiles:",40),export_phi_gen_freq
        write(6  ,'(3X,A40,I9)')adjl("Export density profiles:",40),export_phi_gen_freq
    endif
else
    export_phi_gen_freq = 1
    write(iow,'(3X,A40,I9)')adjl("Export density profiles:",40),export_phi_gen_freq
    write(6  ,'(3X,A40,I9)')adjl("Export density profiles:",40),export_phi_gen_freq
endif


if (log_export_phi_indiv_freq) then
    if (export_phi_indiv_freq.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export individual density profiles:",40),export_phi_indiv_freq
        write(6  ,'(3X,A40,I9)')adjl("Export individual density profiles:",40),export_phi_indiv_freq
    else
        export_phi_indiv_freq = 0
        write(iow,'(3X,A40,I9)')adjl("Export individual density profiles:",40),export_phi_indiv_freq
        write(6  ,'(3X,A40,I9)')adjl("Export individual density profiles:",40),export_phi_indiv_freq
    endif
else
    export_phi_indiv_freq = 0
    write(iow,'(3X,A40,I9)')adjl("Export individual density profiles:",40),export_phi_indiv_freq
    write(6  ,'(3X,A40,I9)')adjl("Export individual density profiles:",40),export_phi_indiv_freq
endif


if (log_export_field_freq) then
    if (export_field_freq.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export field:",40),export_field_freq
        write(6  ,'(3X,A40,I9)')adjl("Export field:",40),export_field_freq
    else
        export_field_freq = 0
        write(iow,'(3X,A40,I9)')adjl("Export field:",40),export_field_freq
        write(6  ,'(3X,A40,I9)')adjl("Export field:",40),export_field_freq
    endif
else
    export_field_freq = 1
    write(iow,'(3X,A40,I9)')adjl("Export field:",40),export_field_freq
    write(6  ,'(3X,A40,I9)')adjl("Export field:",40),export_field_freq
endif


if (log_export_field_bin_freq) then
    if (export_field_bin_freq.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export binary field:",40),export_field_bin_freq
        write(6  ,'(3X,A40,I9)')adjl("Export binary field:",40),export_field_bin_freq
    else
        export_field_bin_freq = 0
        write(iow,'(3X,A40,I9)')adjl("Export binary field:",40),export_field_bin_freq
        write(6  ,'(3X,A40,I9)')adjl("Export binary field:",40),export_field_bin_freq
    endif
else
    export_field_bin_freq = 1
    write(iow,'(3X,A40,I9)')adjl("Export binary field:",40),export_field_bin_freq
    write(6  ,'(3X,A40,I9)')adjl("Export binary field:",40),export_field_bin_freq
endif


if (log_export_propagators_freq) then
    if (export_propagators_freq.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export propagators:",40),export_propagators_freq
        write(6  ,'(3X,A40,I9)')adjl("Export propagators:",40),export_propagators_freq
    else
        export_propagators_freq = 0
        write(iow,'(3X,A40,I9)')adjl("Export propagators:",40),export_propagators_freq
        write(6  ,'(3X,A40,I9)')adjl("Export propagators:",40),export_propagators_freq
    endif
else
    export_propagators_freq = 0
    write(iow,'(3X,A40,I9)')adjl("Export propagators:",40),export_propagators_freq
    write(6  ,'(3X,A40,I9)')adjl("Export propagators:",40),export_propagators_freq
endif


if (log_export_brush_thickness_freq) then
    if (export_brush_thickness_freq.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export brush thickness:",40),export_brush_thickness_freq
        write(6  ,'(3X,A40,I9)')adjl("Export brush thickness:",40),export_brush_thickness_freq
    else
        export_brush_thickness_freq = 0
        write(iow,'(3X,A40,I9)')adjl("Export brush thickness:",40),export_brush_thickness_freq
        write(6  ,'(3X,A40,I9)')adjl("Export brush thickness:",40),export_brush_thickness_freq
    endif
else
    export_brush_thickness_freq = 1
    write(iow,'(3X,A40,I9)')adjl("Export brush thickness:",40),export_brush_thickness_freq
    write(6  ,'(3X,A40,I9)')adjl("Export brush thickness:",40),export_brush_thickness_freq
endif


if (log_export_chains_per_area_freq) then
    if (export_chains_per_area_freq.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export chains per area profiles:",40),export_chains_per_area_freq
        write(6  ,'(3X,A40,I9)')adjl("Export chains per area profiles:",40),export_chains_per_area_freq
    else
        export_chains_per_area_freq = 0
        write(iow,'(3X,A40,I9)')adjl("Export chains per area profiles:",40),export_chains_per_area_freq
        write(6  ,'(3X,A40,I9)')adjl("Export chains per area profiles:",40),export_chains_per_area_freq
    endif
else
    export_chains_per_area_freq = 0
    write(iow,'(3X,A40,I9)')adjl("Export chains per area profiles:",40),export_chains_per_area_freq
    write(6  ,'(3X,A40,I9)')adjl("Export chains per area profiles:",40),export_chains_per_area_freq
endif


if (log_export_ads_free_freq) then
    if (export_ads_free_freq.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40),export_ads_free_freq
        write(6  ,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40),export_ads_free_freq
    else
        export_ads_free_freq = 0
        write(iow,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40),export_ads_free_freq
        write(6  ,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40),export_ads_free_freq
    endif
else
    export_ads_free_freq = 0
    write(iow,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40),export_ads_free_freq
    write(6  ,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40),export_ads_free_freq
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('----------------------------------BOUNDARY CONDITIONS---------------------------------',85)
write(*  ,'(A85)')adjl('----------------------------------BOUNDARY CONDITIONS---------------------------------',85)


if (log_n_dirichlet_faces.and.n_dirichlet_faces>0) then
    write(iow,'(3X,A40,I9)')adjl("Number of Dirichlet (q=0) faces:",40), n_dirichlet_faces
    write(iow,'(3X,A39)',advance='no')"Face ids:                               "
    write(6  ,'(3X,A40,I9)')adjl("Number of Dirichlet (q=0) faces:",40), n_dirichlet_faces
    write(6  ,'(3X,A39)',advance='no')"Face ids:                               "
    do ii = 1, n_dirichlet_faces
        write(iow,'(I3)',advance='no') ids_dirichlet_faces(ii)
        write(6  ,'(I3)',advance='no') ids_dirichlet_faces(ii)
    enddo
    write(iow,*)
    write(6  ,*)
else
    n_dirichlet_faces = 0
    allocate(ids_dirichlet_faces(1))
    ids_dirichlet_faces(1)=-1
    write(iow,'(3X,A40)')adjl("There are no Dirichlet (q=0) faces.",40)
    write(6  ,'(3X,A40)')adjl("There are no Dirichlet (q=0) faces.",40)
endif


if (log_n_nanopart_faces.and.n_nanopart_faces>0) then
    write(iow,'(3X,A40,I9)')adjl("Number of Nanoparticle (q=0) faces:",40), n_nanopart_faces
    write(iow,'(3X,A39)',advance='no')"Face ids:                               "
    write(6  ,'(3X,A40,I9)')adjl("Number of Nanoparticle (q=0) faces:",40), n_nanopart_faces
    write(6  ,'(3X,A39)',advance='no')"Face ids:                               "
    do ii = 1, n_nanopart_faces
        write(iow,'(I3)',advance='no') ids_nanopart_faces(ii)
        write(6  ,'(I3)',advance='no') ids_nanopart_faces(ii)
    enddo
    write(iow,*)
    write(6  ,*)
else
    n_nanopart_faces = 0
    allocate(ids_nanopart_faces(1))
    ids_nanopart_faces(1)=-1
    write(iow,'(3X,A38)')adjl("There are no nanoparticle (q=0) faces.",40)
    write(6  ,'(3X,A38)')adjl("There are no nanoparticle (q=0) faces.",40)
endif


if (log_periodicity.and.periodicity>0) then
    domain_is_periodic = .true.

    if (periodic_axis_id(1)) then
        write(iow,'(3X,A40)')adjl("The domain is periodic along the x-axis",40)
        write(iow,'(6X,A10,2I3)') "Face ids: ", periodic_face_id(1), periodic_face_id(2)
        write(6  ,'(3X,A40)')adjl("The domain is periodic along the x-axis",40)
        write(6  ,'(6X,A10,2I3)') "Face ids: ", periodic_face_id(1), periodic_face_id(2)
    endif
    if (periodic_axis_id(2)) then
        write(iow,'(3X,A40)')adjl("The domain is periodic along the y-axis",40)
        write(iow,'(6X,A10,2I3)') "Face ids: ", periodic_face_id(3), periodic_face_id(4)
        write(6  ,'(3X,A40)')adjl("The domain is periodic along the y-axis",40)
        write(6  ,'(6X,A10,2I3)') "Face ids: ", periodic_face_id(3), periodic_face_id(4)
    endif
    if (periodic_axis_id(3)) then
        write(iow,'(3X,A40)')adjl("The domain is periodic along the z-axis",40)
        write(iow,'(6X,A10,2I3)') "Face ids: ", periodic_face_id(5), periodic_face_id(6)
        write(6  ,'(3X,A40)')adjl("The domain is periodic along the z-axis",40)
        write(6  ,'(6X,A10,2I3)') "Face ids: ", periodic_face_id(5), periodic_face_id(6)
    endif
else
    domain_is_periodic = .false.
    write(iow,'(3X,A38)')adjl("The domain is not periodic.",40)
    write(6  ,'(3X,A38)')adjl("The domain is not periodic.",40)
endif



write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('---------------------------------HAMAKER INTERACTIONS---------------------------------',85)
write(*  ,'(A85)')adjl('---------------------------------HAMAKER INTERACTIONS---------------------------------',85)


if (log_sigma_polymer) then
    if (sigma_pol>0) then
        write(iow,'(3X,A40,E16.9,A11)')adjl("Sigma of polymer:",40),sigma_pol," [Angstrom]"
        write(6  ,'(3X,A40,E16.9,A11)')adjl("Sigma of polymer:",40),sigma_pol," [Angstrom]"
    else
        write(ERROR_MESSAGE,'("sigma_polymer is negative: ",E16.9," Angstroms")') sigma_pol
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="sigma_polymer not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif


if (log_Hamaker_constant_of_polymer) then
    write(iow,'(3X,A40,E16.9,A10)')adjl("Hamaker constant of polymer:",40),A_pol," [10-20 J]"
    write(6  ,'(3X,A40,E16.9,A10)')adjl("Hamaker constant of polymer:",40),A_pol," [10-20 J]"
    A_pol = A_pol * 1e-20
else
    A_pol = 0.d0
    write(iow,'(3X,A40)')adjl("Hamaker constant for pol not found.",40)
    write(iow,'(3X,A40,E16.9,A10)')adjl("It was set to the default value: ",40),A_pol," [10-20 J]"
    write(6  ,'(3X,A40)')adjl("Hamaker constant for pol not found.",40)
    write(6  ,'(3X,A40,E16.9,A10)')adjl("It was set to the default value: ",40),A_pol," [10-20 J]"
endif


if (log_wall_distance) then
    write(iow,'(3X,A40,E16.9,A11)')adjl("Wall distance for Hamaker:",40),wall_distance," [Angstrom]"
    write(6  ,'(3X,A40,E16.9,A11)')adjl("Wall distance for Hamaker:",40),wall_distance," [Angstrom]"
else
    wall_distance = 5.d0
    write(iow,'(3X,A24)')adjl("Wall distance not found.",40)
    write(iow,'(3X,A40,E16.9,A11)')adjl("It was set to the default value: ",40),wall_distance," [Angstrom]"
    write(6  ,'(3X,A24)')adjl("Wall distance not found.",40)
    write(6  ,'(3X,A40,E16.9,A11)')adjl("It was set to the default value: ",40),wall_distance," [Angstrom]"
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-------------------------------NONBONDED INTERACTIONS---------------------------------',85)
write(*  ,'(A85)')adjl('-------------------------------NONBONDED INTERACTIONS---------------------------------',85)


if (log_eos_type) then
    if (eos_type.eq.eos_helfand) then
        write(iow,'(3X,A40,I9)')adjl("Equation of state: Helfand",40),eos_type
        write(*  ,'(3X,A40,I9)')adjl("Equation of state: Helfand",40),eos_type
    elseif (eos_type.eq.eos_sl) then
        write(iow,'(3X,A45,I9)')adjl("Equation of state: Sanchez-Lacombe",40),eos_type
        write(*  ,'(3X,A45,I9)')adjl("Equation of state: Sanchez-Lacombe",40),eos_type
    else
        write(*  ,'(A45,I11)') 'EOS flag different than 0 (HF) or 1 (SL)',eos_type
        STOP
    endif
    call init_scf_params()
else
    write(iow,'(A45)') 'EOS flag not set'
    write(*  ,'(A45)') 'EOS flag not set'
    STOP
endif


if (log_eos_coeffs) then
    if (eos_type.eq.eos_helfand) then
        write(iow,'(3X,A40,E16.9,A8)')adjl("Helfand isothermal compressibility:",40),hlf_kappa_T," [Pa^-1]"
        write(*  ,'(3X,A40,E16.9,A8)')adjl("Helfand isothermal compressibility:",40),hlf_kappa_T," [Pa^-1]"
    elseif (eos_type.eq.eos_sl) then
        write(iow,'(A40,3(F16.4))') "rho_star, T_star, P_star = ", rho_star, T_star, P_star
        write(*  ,'(A40,3(F16.4))') "rho_star, T_star, P_star = ", rho_star, T_star, P_star
    endif
else
    write(iow,'(A40)') "EOS coeffs were not found"
    write(*  ,'(A40)') "EOS coeffs were not found"
    STOP
endif

if  (log_influence_param) then
    square_gradient = .true.
    write(iow,'(3X,A40,E16.9,A14)')adjl("Influence parameter:",45), k_gr_tilde, " [J*m^5/mol^2]"
    write(*  ,'(3X,A40,E16.9,A14)')adjl("Influence parameter:",45), k_gr_tilde, " [J*m^5/mol^2]"
else
    k_gr_tilde = 0.d0
    square_gradient = .false.
    write(iow,'(3X,A40,E16.9,A14)')adjl("Influence parameter not found. Auto:",45), k_gr_tilde, " [J*m^5/mol^2]"
    write(*  ,'(3X,A40,E16.9,A14)')adjl("Influence parameter not found. Auto:",45), k_gr_tilde, " [J*m^5/mol^2]"
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-------------------------------------SPATIAL MESH------------------------------------',85)
write(*  ,'(A85)')adjl('-------------------------------------SPATIAL MESH------------------------------------',85)


if (log_mesh_filename) then
    inquire(file=mesh_filename, exist=file_exists)
    if (.not.file_exists) then
        write(ERROR_MESSAGE,'("Mesh file ",A16," does not exist!")')mesh_filename
        call exit_with_error(1,1,1,ERROR_MESSAGE)
        STOP
    endif
    write(iow,'(3X,A40,A16)')adjl("Reading mesh from file:",40),mesh_filename
    write(6  ,'(3X,A40,A16)')adjl("Reading mesh from file:",40),mesh_filename
else
    mesh_filename = "./in.mesh"
    write(iow,'(3X,A40)')adjl("Mesh input file not specified.",40)
    write(iow,'(3X,A40,A16)')adjl("Reading default mesh file:",40),mesh_filename
    write(6  ,'(3X,A40)')adjl("Mesh input file not specified.",40)
    write(6  ,'(3X,A40,A16)')adjl("Reading default mesh file:",40),mesh_filename

    inquire(file=mesh_filename, exist=file_exists)
    if (.not.file_exists) then
        write(ERROR_MESSAGE,'("Default mesh file ",A16," does not exist!")')mesh_filename
        call exit_with_error(1,1,1,ERROR_MESSAGE)
        STOP
    endif
endif

return
!--------------------------------------------------------------------------------!
end subroutine parser_input
