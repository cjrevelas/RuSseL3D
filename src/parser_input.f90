!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine parser_input()
!--------------------------------------------------------------------------------!
use parser_vars_mod
use write_helper_mod
use error_handing_mod
use flags_mod
use eos_mod
use iofiles_mod
use defaults_mod
!--------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------!
character(100) :: line
character(1)   :: axis

integer :: Reason, ii, id, face1, face2

real(8) :: dirichletValue = dflt_dirichletValue
real(8) :: nanoparticleValue = dflt_nanoparticleValue

logical :: log_meshFile                    = .False.
logical :: log_graftFile                   = .False.
logical :: log_fieldFile                   = .False.
logical :: log_temperature                 = .False.
logical :: log_pressure                    = .False.
logical :: log_massDensity                 = .False.
logical :: log_monomer_mass                = .False.
logical :: log_number_of_iterations        = .False.
logical :: log_set_initial_iteration       = .False.
logical :: log_fieldTol                    = .False.
logical :: log_freeEnergyTol               = .False.
logical :: log_freeEnergyTolForDelta       = .False.
logical :: log_field_init_scheme           = .False.
logical :: log_fraction_of_new_field       = .False.
logical :: log_mx_exist                    = .False.
logical :: log_rg2OfMatrixMonomer          = .False.
logical :: log_ds_ave_mx                   = .False.
logical :: log_lengthMatrix                = .False.
logical :: log_contour_discr_mx            = .False.
logical :: log_ads_distance                = .False.
logical :: log_gr_exist                    = .False.
logical :: log_rg2OfGraftedMonomer         = .False.
logical :: log_lengthGrafted               = .False.
logical :: log_ds_ave_gr                   = .False.
logical :: log_contour_discr_gr            = .False.
logical :: log_grafted_ic_from_delta       = .False.
logical :: log_numGraftedChainsTol         = .False.
logical :: log_r_gpoint                    = .False.
logical :: log_calc_delta_every            = .False.
logical :: log_numDirichletFaces           = .False.
logical :: log_numNanoparticleFaces        = .False.
logical :: log_periodicity                 = .False.
logical :: log_sigma_polymer               = .False.
logical :: log_Hamaker_constant_of_polymer = .False.
logical :: log_wall_distance               = .False.
logical :: log_eos_type                    = .False.
logical :: log_eos_coeffs                  = .False.
logical :: log_influence_param             = .False.
logical :: log_profile_dimension           = .False.
logical :: log_mumpsMatrixType             = .False.
logical :: log_binThickness                = .False.
logical :: log_exportPhiGeneral            = .False.
logical :: log_exportPhiIndividual         = .False.
logical :: log_exportField                 = .False.
logical :: log_exportFieldBinary           = .False.
logical :: log_exportPropagators           = .False.
logical :: log_exportBrushThickness        = .False.
logical :: log_exportChainsPerArea         = .False.
logical :: log_exportAdsorbedFree          = .False.
!--------------------------------------------------------------------------------!
! Initialize periodicity arrays
do ii = 1, 3
    periodicAxisId = .False.
enddo
do ii = 1, 6
    periodicFaceId = -1
enddo

! Parse input file to retrieve simulation parameters
call GET_COMMAND_ARGUMENT(1,inputFile)

if (inputFile == '') inputFile = dflt_inputFile
inquire(file=inputFile, exist=file_exists)

if (file_exists) then
    open(unit=256, file = inputFile)
else
    write(ERROR_MESSAGE,'("File ",A15," does not exist!")') inputFile
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
            read(line,'(A40)') meshFile
            log_meshFile = .True.
        elseif (INDEX(line,"# gpoints input file") > 0) then
            read(line,'(A40)') graftFile
            log_graftFile = .True.
        elseif (INDEX(line,"# field input file") > 0) then
            read(line,*) fieldFile
            log_fieldFile = .True.
        elseif (INDEX(line,"# temp") > 0) then
            read(line,*) temperature
            log_temperature = .True.
        elseif (INDEX(line," pres") > 0) then
            read(line,*) pressure
            log_pressure = .True.
        elseif (INDEX(line,"# mass den") > 0) then
            read(line,*) massDensity
            log_massDensity = .True.
        elseif (INDEX(line,"# mon mass") > 0) then
            read(line,*) massOfMonomer
            log_monomer_mass = .True.
        elseif (INDEX(line,"# Rg2/mon matrix") > 0) then
            read(line,*) rg2OfMatrixMonomer
            log_rg2OfMatrixMonomer = .True.
        elseif (INDEX(line,"# Rg2/mon grafted") > 0) then
            read(line,*) rg2OfGraftedMonomer
            log_rg2OfGraftedMonomer = .True.
        elseif (INDEX(line,"# chain length matrix") > 0) then
            read(line,*) lengthMatrix
            log_lengthMatrix = .True.
        elseif (INDEX(line,"# chain length grafted") > 0) then
            read(line,*) lengthGrafted
            log_lengthGrafted = .True.
        elseif (INDEX(line,"# num gr chains tol") > 0) then
            read(line,*) numGraftedChainsTol
            log_numGraftedChainsTol = .True.
        elseif (INDEX(line,"# pol sigma") > 0) then
            read(line,*) sigma_pol
            log_sigma_polymer = .True.
        elseif (INDEX(line,"# pol ham") > 0) then
            read(line,*) A_pol
            log_Hamaker_constant_of_polymer = .True.
        elseif (INDEX(line,"# wall dist") > 0) then
            read(line,*) wall_distance
            log_wall_distance = .True.
        elseif (INDEX(line,"# fraction") > 0) then
            read(line,*) frac
            log_fraction_of_new_field = .True.
        elseif (INDEX(line,"# max error") > 0) then
            read(line,*) fieldTol
            log_fieldTol = .True.
        elseif (INDEX(line,"# max energy error") > 0) then
            read(line,*) freeEnergyTol
            log_freeEnergyTol = .True.
        elseif (INDEX(line,"# energy tol for delta") > 0) then
            read(line,*) freeEnergyTolForDelta
            log_freeEnergyTolForDelta = .True.
        elseif (INDEX(line,"# init iter") > 0) then
            read(line,*) init_iter
            log_set_initial_iteration= .True.
        elseif (INDEX(line,"# num iter") > 0) then
            read(line,*) iterations
            log_number_of_iterations = .True.
        elseif (INDEX(line,"# init field") > 0) then
            read(line,*) field_init_scheme
            log_field_init_scheme = .True.
        elseif (INDEX(line,"# bin thickness") > 0) then
            read(line,*) binThickness
            log_binThickness = .True.
        elseif (INDEX(line,"# export dens profs") > 0) then
            read(line,*) exportPhiGeneral
            log_exportPhiGeneral = .True.
        elseif (INDEX(line,"# export indiv dens profs") > 0) then
            read(line,*) exportPhiIndividual
            log_exportPhiIndividual = .True.
        elseif (INDEX(line,"# export field") > 0) then
            read(line,*) exportField
            log_exportField = .True.
        elseif (INDEX(line,"# export binary field") > 0) then
            read(line,*) exportFieldBinary
            log_exportFieldBinary = .True.
        elseif (INDEX(line,"# export propagators") > 0) then
            read(line,*) exportPropagators
            log_exportPropagators = .True.
        elseif (INDEX(line,"# export brush thickness") > 0) then
            read(line,*) exportBrushThickness
            log_exportBrushThickness = .True.
        elseif (INDEX(line,"# export chains per area profs") > 0) then
            read(line,*) exportChainsPerArea
            log_exportChainsPerArea = .True.
        elseif (INDEX(line,"# export ads vs free profs") > 0) then
            read(line,*) exportAdsorbedFree
            log_exportAdsorbedFree = .True.
        elseif (INDEX(line,"# prof dim") > 0) then
            read(line,*) prof_dim
            log_profile_dimension = .True.
        elseif (INDEX(line,"# use matrix") > 0) then
            read(line,*) mx_exist
            log_mx_exist = .True.
        elseif (INDEX(line,"# use grafted") > 0) then
            read(line,*) gr_exist
            log_gr_exist = .True.
        elseif (INDEX(line,"# contour step matrix") > 0) then
            read(line,*) ds_ave_mx_ed, ds_ave_mx_conv
            log_ds_ave_mx = .True.
        elseif (INDEX(line,"# discr scheme matrix") > 0) then
            read(line,*) contour_discr_mx
            if (contour_discr_mx.eq.contour_hybrid) read(256,*) xs_crit_mx
            log_contour_discr_mx = .True.
        elseif (INDEX(line,"# ads distance") > 0) then
            read(line,*) ads_distance
            log_ads_distance = .True.
        elseif (INDEX(line,"# contour step grafted") > 0) then
            read(line,*) ds_ave_gr_ed, ds_ave_gr_conv
            log_ds_ave_gr = .True.
        elseif (INDEX(line,"# gp dist from solid") > 0) then
            read(line,*) r_gpoint
            log_r_gpoint = .True.
        elseif (INDEX(line,"# discr scheme grafted") > 0) then
            read(line,*) contour_discr_gr
            if (contour_discr_gr.eq.contour_hybrid) read(256,*) xs_crit_gr
            log_contour_discr_gr = .True.
        elseif (INDEX(line,"# mumps matrix") > 0) then
            read(line,*) mumpsMatrixType
            log_mumpsMatrixType = .True.
        elseif (INDEX(line,"# eos type") > 0) then
            read(line,*) eos_type
            log_eos_type = .True.
        elseif (INDEX(line,"# eos coeffs") > 0) then
            if (eos_type.eq.eos_helfand) then
                read(line,*) hlf_kappa_T
            elseif (eos_type.eq.eos_sl)  then
                read(line,*) rho_star, T_star, P_star
            endif
            log_eos_coeffs = .True.
        elseif (INDEX(line, "# eos infl param") > 0) then
            read(line,*) k_gr_tilde
            log_influence_param = .True.
        elseif (INDEX(line,"# calc delta") > 0) then
            read(line,*) grafted_ic_from_delta
            log_grafted_ic_from_delta = .True.
        elseif (INDEX(line,"# freq delta") > 0) then
            read(line,*) calc_delta_every
            log_calc_delta_every = .True.
        elseif (INDEX(line,"# num faces") > 0) then
            read(line,*) numDirichletFaces
            if (numDirichletFaces > 0) then
                allocate(dirichletFaceId(numDirichletFaces))
                allocate(dirichletFaceValue(numDirichletFaces))
                allocate(sigma_plate(numDirichletFaces))
                allocate(A_plate(numDirichletFaces))
                do ii = 1, numDirichletFaces
                    read(256,*) id, sigma_plate(ii), A_plate(ii), dirichletValue
                    dirichletFaceId(ii) = id
                    dirichletFaceValue(ii) = dirichletValue
                enddo
                A_plate = A_plate * 1.0d-20
                log_numDirichletFaces = .True.
            endif
        elseif (INDEX(line,"# num nanop") > 0) then
            read(line,*) numNanoparticleFaces
            if (numNanoparticleFaces > 0) then
                allocate(nanoparticleFaceId(numNanoparticleFaces))
                allocate(nanoparticleFaceValue(numNanoparticleFaces))
                allocate(center_np(3,numNanoparticleFaces))
                allocate(radius_np_eff(numNanoparticleFaces))
                allocate(sigma_np(numNanoparticleFaces))
                allocate(A_np(numNanoparticleFaces))
                do ii = 1, numNanoparticleFaces
                    read(256,*) id, radius_np_eff(ii), center_np(1,ii), center_np(2,ii), center_np(3,ii), &
                                                          & sigma_np(ii), A_np(ii), nanoparticleValue
                    nanoparticleFaceId(ii) = id
                    nanoparticleFaceValue(ii) = nanoparticleValue
                enddo
                A_np = A_np * 1.0d-20
                log_numNanoparticleFaces = .True.
            endif
        elseif (INDEX(line,"# periodicity") > 0) then
            read(line,*) periodicity
            do ii = 1, periodicity
                read(256,*) axis, face1, face2
                if (axis=='x') then
                    periodicAxisId(1) = .True.
                    periodicFaceId(1) = face1
                    periodicFaceId(2) = face2
                elseif (axis=='y') then
                    periodicAxisId(2) = .True.
                    periodicFaceId(3) = face1
                    periodicFaceId(4) = face2
                elseif (axis=='z') then
                    periodicAxisId(3) = .True.
                    periodicFaceId(5) = face1
                    periodicFaceId(6) = face2
                else
                    write(ERROR_MESSAGE,'("Periodicity axis is not valid. Must be x, y or z.")')
                    call exit_with_error(1,1,1,ERROR_MESSAGE)
                endif
            enddo
            log_periodicity = .True.
        endif
   endif
enddo

close(256)

! Check input parameters
write(iow,'(A85)')adjl('-----------------------------------SYSTEM PARAMETERS---------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------SYSTEM PARAMETERS---------------------------------',85)


if (log_temperature) then
    if (temperature>0.0d0) then
        write(iow,'(3X,A40,E16.9,A4)')adjl("Temperature:",40), temperature," [K]"
        write(6  ,'(3X,A40,E16.9,A4)')adjl("Temperature:",40), temperature," [K]"
        beta = 1.0d0 / (boltz_const_Joule_K * temperature)
    else
        write(ERROR_MESSAGE,'("Temperature is negative: ",E16.9," K")') temperature
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE = "Temperature was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif


if (log_pressure) then
    if (pressure>=0.0d0) then
        write(iow,'(3X,A40,E16.9,A6)')adjl("Pressure:",40), pressure, " [atm]"
        write(*  ,'(3X,A40,E16.9,A6)')adjl("Pressure:",40), pressure, " [atm]"
        pressure = pressure * atm_to_pa
    else
        write(ERROR_MESSAGE,'("Pressure is negative: ",E16.9, " atm")') pressure
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    pressure = dflt_pressure
    write(iow,'(A40)') "Pressure was set to 0 atm"
    write(*  ,'(A40)') "Pressure was set to 0 atm"
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('----------------------------------POLYMER PROPERTIES---------------------------------',85)
write(*  ,'(A85)')adjl('----------------------------------POLYMER PROPERTIES---------------------------------',85)


if (log_massDensity) then
    if (massDensity>0.0d0) then
        write(iow,'(3X,A40,E16.9,A8)')adjl("Mass density:",40), massDensity, " [g/cm3]"
        write(6  ,'(3X,A40,E16.9,A8)')adjl("Mass density:",40), massDensity, " [g/cm3]"
    else
        write(ERROR_MESSAGE,'("Mass density is negative: ",E16.9," g/cm3")') massDensity
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE = "Mass density was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif


if (log_monomer_mass) then
    if (massOfMonomer>0.0d0) then
        write(iow,'(3X,A40,E16.9,A8)')adjl("Monomer mass:",40), massOfMonomer, "[g/mol]"
        write(6  ,'(3X,A40,E16.9,A8)')adjl("Monomer mass:",40), massOfMonomer, "[g/mol]"
    else
        write(ERROR_MESSAGE,'("Monomer mass is negative: ",E16.9)') massOfMonomer
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE = "Monomer mass not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif


if (mx_exist.eq.1) then
    if (log_lengthMatrix) then
        if (lengthMatrix>0.0d0) then
            lengthMatrixMax = lengthMatrix
        else
            write(ERROR_MESSAGE,'("Chain length of matrix chains is negative: ",E16.9)') lengthMatrix
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE = "Chain length of matrix chains was not detected."
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
    if (log_graftFile) then
        inquire(file=graftFile, exist=file_exists)
        if (.not.file_exists) then
            write(ERROR_MESSAGE,'("Grafting points file ",A16," does not exist!")') graftFile
            call exit_with_error(1,1,1,ERROR_MESSAGE)
            STOP
        endif
        write(iow,'(3X,A40,A16)')adjl("Reading grafting points from file:",40), graftFile
        write(6  ,'(3X,A40,A16)')adjl("Reading grafting points from file:",40), graftFile
    else
        graftFile = dflt_graftFile
        write(iow,'(3X,A40)')adjl("Grafting points input file not specified.",40)
        write(iow,'(3X,A40,A16)')adjl("Reading default grafting points file:",40), graftFile
        write(6  ,'(3X,A40)')adjl("Grafting points input file not specified.",40)
        write(6  ,'(3X,A40,A16)')adjl("Reading default grafting points file:",40), graftFile


        inquire(file=graftFile, exist=file_exists)
        if (.not.file_exists) then
            write(ERROR_MESSAGE,'("Default grafting points file ",A16," does not exist!")') graftFile
            call exit_with_error(1,1,1,ERROR_MESSAGE)
            STOP
        endif
    endif

    if (log_rg2OfGraftedMonomer) then
        if (rg2OfGraftedMonomer>0.0d0) then
            write(iow,'(3X,A40,E16.9,A13)')adjl("Rg2 per grafted monomer:",40), rg2OfGraftedMonomer, "[Angstrom^2]"
            write(6  ,'(3X,A40,E16.9,A13)')adjl("Rg2 per grafted monomer:",40), rg2OfGraftedMonomer, "[Angstrom^2]"
        else
            write(ERROR_MESSAGE,'("Rg2 per grafted monomer is negative: ",E16.9)') rg2OfGraftedMonomer
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE = "Rg2 per matrix monomer was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_lengthGrafted) then
        if (lengthGrafted>0.0d0) then
            write(iow,'(3X,A40,E16.9,A11)')adjl("Chain length of grafted chains:",40), lengthGrafted,"[monomers]"
            write(iow,'(3X,A40,E16.9,A11)')adjl("Radius of gyration of grafted chains:",40), &
                                                              & DSQRT(rg2OfGraftedMonomer*lengthGrafted), "[Angstrom]"
            write(6  ,'(3X,A40,E16.9,A11)')adjl("Chain length of grafted chains:",40), lengthGrafted, "[monomers]"
            write(6  ,'(3X,A40,E16.9,A11)')adjl("Radius of gyration of grafted chains:",40), &
                                                              & DSQRT(rg2OfGraftedMonomer*lengthGrafted), "[Angstrom]"
            lengthMatrixMax = MAX(lengthMatrix, lengthGrafted)
        else
            write(ERROR_MESSAGE,'("Chain length of grafted chains is negative: ",E16.9)') lengthGrafted
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE = "Chain length of grafted chains was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_ds_ave_gr) then
        if (ds_ave_gr_ed>0 .and. ds_ave_gr_conv>0) then
            ns_gr_ed   = 2 * NINT(0.5d0 * lengthGrafted / ds_ave_gr_ed)
            ns_gr_conv = 2 * NINT(0.5d0 * lengthGrafted / ds_ave_gr_conv)
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
        ERROR_MESSAGE = "Contour step of grafted chains was not detected."
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
                ERROR_MESSAGE = "Critical contour point of grafted chains was not detected."
                call exit_with_error(1,1,1,ERROR_MESSAGE)
            endif

            call compute_contour_step(ds_ave_gr_ed, xs_crit_gr, lengthGrafted, ns_gr_ed)

        elseif (contour_discr_gr.eq.contour_asymm) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "asymm"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of graftd chains:", "asymm"
        else
            write(ERROR_MESSAGE,'("Not valid Edwards contour scheme of grafted chains: ",I5)') contour_discr_gr
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE = "Edwards contour scheme of grafted chains was not detected."
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
    if (log_rg2OfMatrixMonomer) then
        if (rg2OfMatrixMonomer>0.0d0) then
            write(iow,'(3X,A40,E16.9,A13)')adjl("Rg2 per matrix monomer:",40), rg2OfMatrixMonomer, "[Angstrom^2]"
            write(6  ,'(3X,A40,E16.9,A13)')adjl("Rg2 per matrix monomer:",40), rg2OfMatrixMonomer, "[Angstrom^2]"
        else
            write(ERROR_MESSAGE,'("Rg2 per matrix monomer is negative: ",E16.9)') rg2OfMatrixMonomer
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Rg2 per matrix monomer was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    write(iow,'(3X,A40,E16.9,A11)')adjl("Chain length of matrix chains:",40), lengthMatrix, "[monomers]"
    write(iow,'(3X,A40,E16.9,A11)')adjl("Radius of gyration of matrix chains:",40), &
                                                             & DSQRT(rg2OfMatrixMonomer*lengthMatrix), "[Angstrom]"
    write(6  ,'(3X,A40,E16.9,A11)')adjl("Chain length of matrix chains:",40), lengthMatrix, "[monomers]"
    write(6  ,'(3X,A40,E16.9,A11)')adjl("Radius of gyration of matrix chains:",40), &
                                                             & DSQRT(rg2OfMatrixMonomer*lengthMatrix), "[Angstrom]"

    if (log_ds_ave_mx) then
        if (ds_ave_mx_ed>0.0d0 .and. ds_ave_mx_conv>0.0d0) then
            ns_mx_ed   = 2 * NINT(0.5d0 * lengthMatrixMax / ds_ave_mx_ed)
            ns_mx_conv = 2 * NINT(0.5d0 * lengthMatrix    / ds_ave_mx_conv)

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
                ERROR_MESSAGE = "Critical contour point of matrix chains was not detected."
                call exit_with_error(1,1,1,ERROR_MESSAGE)
            endif

            call compute_contour_step(ds_ave_mx_ed, xs_crit_mx, lengthMatrixMax, ns_mx_ed)

        elseif (contour_discr_mx.eq.contour_asymm) then
            write(iow,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "asymm"
            write(6  ,'(3X,A40,A16)') "Edwards contour scheme of matrix chains:", "asymm"
        else
            write(ERROR_MESSAGE,'("Not valid Edwards contour scheme of matrix chains: ",I5)') contour_discr_mx
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE = "Edwards contour scheme of matrix chains was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_ads_distance) then
        if (ads_distance.ge.0.0d0) then
            write(iow,'(3X,A40,E16.9,A11)')adjl("Adsorption distance for chain segments:",40), ads_distance, "[Angstrom]"
            write(6  ,'(3X,A40,E16.9,A11)')adjl("Adsorption distance for chain segments:",40), ads_distance, "[Angstrom]"
        else
            write(ERROR_MESSAGE,'("Adsorption distance for chain segments is negative:",E16.9)') ads_distance
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ads_distance = dflt_ads_distance
        write(iow,'(3X,A40)')adjl("Adsorption distance not found.",40)
        write(iow,'(3X,A40,E16.9,A11)')adjl("It was set to the default value:",40), ads_distance, "[Angstrom]"
        write(6  ,'(3X,A40)')adjl("Adsorption distance not found.",40)
        write(6  ,'(3X,A40,E16.9,A11)')adjl("It was set to the default value:",40), ads_distance, "[Angstrom]"
    endif
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('--------------------------------SIMULATION PARAMETERS---------------------------------',85)
write(*  ,'(A85)')adjl('--------------------------------SIMULATION PARAMETERS---------------------------------',85)


if (log_mx_exist.and.mx_exist.ge.1) then
    mx_exist = 1
    write(iow,'(3X,A40,I9)')adjl("The system includes matrix chains:",40), mx_exist
    write(6  ,'(3X,A40,I9)')adjl("The system includes matrix chains:",40), mx_exist
else
    mx_exist = dflt_mx_exist
    write(iow,'(3X,A40,I9)')adjl("System does not include matrix chains:",40), mx_exist
    write(6  ,'(3X,A40,I9)')adjl("System does not include matrix chains:",40), mx_exist
endif


if (log_gr_exist.and.gr_exist.ge.1) then
    gr_exist = 1
    write(iow,'(3X,A40,I9)')adjl("The system includes grafted chains:",40), gr_exist
    write(6  ,'(3X,A40,I9)')adjl("The system includes grafted chains:",40), gr_exist

    if (log_grafted_ic_from_delta) then
        write(iow,'(3X,A40,I9)')adjl("Grafted ic from delta:",40), grafted_ic_from_delta
        write(6  ,'(3X,A40,I9)')adjl("Grafted ic from delta:",40), grafted_ic_from_delta
        if (grafted_ic_from_delta.ge.1) then
            if (log_calc_delta_every) then
                if (calc_delta_every.gt.0) then
                    write(iow,'(3X,A40,I9)')adjl("Delta is calculated every:",40), calc_delta_every
                    write(6  ,'(3X,A40,I9)')adjl("Delta is calculated every:",40), calc_delta_every
                else
                    write(iow,'(3X,A40)')adjl("Delta is read from file",40)
                    write(6  ,'(3X,A40)')adjl("Delta is read from file",40)
                endif
            endif
            if (log_numGraftedChainsTol) then
                if (numGraftedChainsTol.ge.0.0d0) then
                    write(iow,'(3X,A40,E16.9)')adjl("Number of grafted chains tolerance:",40), numGraftedChainsTol
                    write(6  ,'(3X,A40,E16.9)')adjl("Number of grafted chains tolerance:",40), numGraftedChainsTol
                else
                    write(ERROR_MESSAGE,'("Number of grafted chains tolerance is negative:",E16.9)') numGraftedChainsTol
                    call exit_with_error(1,1,1,ERROR_MESSAGE)
                endif
            else
                numGraftedChainsTol = dflt_numGraftedChainsTol
                write(iow,'(3X,A40)')adjl("Num grafted chains tolerance not found",40)
                write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), numGraftedChainsTol
                write(6  ,'(3X,A40)')adjl("Num grafted chains tolerance not found",40)
                write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), numGraftedChainsTol
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
        if (r_gpoint.gt.0.0d0) then
            write(iow,'(3X,A40,E16.9)')adjl("Distance of grafting points from solid:",40), r_gpoint
            write(6  ,'(3X,A40,E16.9)')adjl("Distance of grafting points from solid:",40), r_gpoint
        else
            write(ERROR_MESSAGE,'("Distance of grafting points from solid must have a positive value:",E16.9)') r_gpoint
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        r_gpoint = dflt_r_gpoint
        write(iow,'(3X,A40)')adjl("Distance of grafting points not found.",40)
        write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), r_gpoint
        write(6  ,'(3X,A40)')adjl("Distance of grafting points not found.",40)
        write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), r_gpoint
    endif
else
    gr_exist = dflt_gr_exist
    write(iow,'(3X,A40,1x,I15)')adjl("System does not include grafted chains:",40), gr_exist
    write(6  ,'(3X,A40,1x,I15)')adjl("System does not include grafted chains:",40), gr_exist
endif


if (log_set_initial_iteration) then
    if (init_iter.lt.0) then
        write(ERROR_MESSAGE,'("Wrong value of initial iteration.",I16)') init_iter
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    init_iter = dflt_init_iter
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
        write(iow,'(3X,A40,I9)')adjl("Maximum number of iterations:",40), iterations
        write(6  ,'(3X,A40,I9)')adjl("Maximum number of iterations:",40), iterations
    else
        write(ERROR_MESSAGE,'("Maximum number of iterations is negative:",I10)') iterations
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    iterations = dflt_iterations
    write(iow,'(3X,A40)')adjl("Max number of iter not found.",40)
    write(iow,'(3X,A40,I9)')adjl("It was set to the default value:",40), iterations
    write(6  ,'(3X,A40)')adjl("Max number of iter not found.",40)
    write(6  ,'(3X,A40,I9)')adjl("It was set to the default value:",40), iterations
endif


if (log_binThickness) then
    if (binThickness.ge.0.0d0) then
        write(iow,'(3X,A40,E16.9,A11)')adjl("Bin thickness:",40), binThickness, "[Angstrom]"
        write(6  ,'(3X,A40,E16.9,A11)')adjl("Bin thickness:",40), binThickness, "[Angstrom]"
    else
        write(ERROR_MESSAGE,'("Bin thickness is negative:",E16.9)') binThickness
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    binThickness = dflt_binThickness
    write(iow,'(3X,A40)')adjl("Bin thickness not found.",40)
    write(iow,'(3X,A40,E16.9,A11)')adjl("It was set to the default value:",40), binThickness, "[Angstrom]"
    write(6  ,'(3X,A40)')adjl("Bin thickness not found.",40)
    write(6  ,'(3X,A40,E16.9,A11)')adjl("It was set to the default value:",40), binThickness, "[Angstrom]"
endif


if (log_profile_dimension) then
    if ( (prof_dim.ge.1).and.(prof_dim.le.3) ) then
        write(iow,'(3X,A40,I9)')adjl("Profile dimension:",40), prof_dim
        write(6  ,'(3X,A40,I9)')adjl("Profile dimension:",40), prof_dim
    else
        write(ERROR_MESSAGE,'("Profile dimension is not between 1 and 3:",I16)') prof_dim
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    prof_dim = dflt_prof_dim
    write(iow,'(3X,A40)')adjl("Profile dimension not found.",40)
    write(iow,'(3X,A40,I8)')adjl("It was set to the default value:",40), prof_dim
    write(6  ,'(3X,A40)')adjl("Profile dimension not found.",40)
    write(6  ,'(3X,A40,I8)')adjl("It was set to the default value:",40), prof_dim
endif


if (log_field_init_scheme) then
    if (field_init_scheme==0) then
        write(iow,'(3X,A34)')adjl("Field will be initialized to zero.",40)
        write(6  ,'(3X,A34)')adjl("Field will be initialized to zero.",40)
    elseif (field_init_scheme==1) then
        if (log_fieldFile) then
            inquire(file=fieldFile, exist=file_exists)
            if (.not.file_exists) then
                write(ERROR_MESSAGE,'("Field input file ",A16," does not exist!")') fieldFile
                call exit_with_error(1,1,1,ERROR_MESSAGE)
                STOP
            endif
            write(iow,'(A43,A16)')adjl("Field will be read from file:",40), fieldFile
            write(6  ,'(A43,A16)')adjl("Field will be read from file:",40), fieldFile
        else
            fieldFile = dflt_fieldFile
            write(iow,'(3X,A40)')adjl("Field input file not specified.",40)
            write(iow,'(3X,A40,A16)')adjl("Reading default field input file:",40), fieldFile
            write(6  ,'(3X,A40)')adjl("Field input file not specified.",40)
            write(6  ,'(3X,A40,A16)')adjl("Reading default field input file:",40), fieldFile

            inquire(file=fieldFile, exist=file_exists)
            if (.not.file_exists) then
                write(ERROR_MESSAGE,'("Default field input file ",A16," does not exist!")') fieldFile
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
else
    write(iow,'(A40)')adjl("Field will be initialized to zero.",40)
    write(6  ,'(A40)')adjl("Field will be initialized to zero.",40)
endif


if (log_mumpsMatrixType) then
    if (mumpsMatrixType == mumps_asymm) then
        write(iow,'(3X,A40,1X,A14,I1,A1)')adjl("MUMPS matrix type:",40), "nonsymmetric (", mumpsMatrixType, ")"
        write(6  ,'(3X,A40,1X,A14,I1,A1)')adjl("MUMPS matrix type:",40), "nonsymmetric (", mumpsMatrixType, ")"
    elseif (mumpsMatrixType == mumps_posDef) then
        write(iow,'(3X,A40,1X,A21,I1,A1)')adjl("MUMPS matrix type:",40), "symmetric pos. def. (", mumpsMatrixType, ")"
        write(6  ,'(3X,A40,1X,A21,I1,A1)')adjl("MUMPS matrix type:",40), "symmetric pos. def. (", mumpsMatrixType, ")"
    elseif (mumpsMatrixType == mumps_genSymm) then
        write(iow,'(3X,A40,1X,A18,I1,A1)')adjl("MUMPS matrix type:",40), "general symmetric(", mumpsMatrixType, ")"
        write(6  ,'(3X,A40,1X,A18,I1,A1)')adjl("MUMPS matrix type:",40), "general symmetric(", mumpsMatrixType, ")"
    else
        write(ERROR_MESSAGE,'("Incorrect MUMPS matrix type.",I16)') mumpsMatrixType
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    mumpsMatrixType = dflt_mumpsMatrixType
    write(iow,'(3X,A40)')adjl("MUMPS matrix type not detected.",40)
    write(iow,'(3X,A40,I10)')adjl("It was set to the nonsymmetric:",40), mumpsMatrixType
    write(6  ,'(3X,A40)')adjl("MUMPS matrix type not detected.",40)
    write(6  ,'(3X,A40,I10)')adjl("It was set to the nonsymmetric:",40), mumpsMatrixType
endif


if (log_fraction_of_new_field) then
    if (frac.ge.0.0d0 .and. frac.le.1.0d0) then
        write(iow,'(3X,A40,E16.9)')adjl("Initial fraction of new field:",40), frac
        write(6  ,'(3X,A40,E16.9)')adjl("Initial fraction of new field:",40), frac
    else
        write(ERROR_MESSAGE,'("Initial fraction of new field is negative or larger than unity:",E16.9)') frac
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    frac = dflt_fraction
    write(iow,'(3X,A40)')adjl("No initial fraction of new field.",40)
    write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), frac
    write(6  ,'(3X,A40)')adjl("No initial fraction of new field.",40)
    write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), frac
endif


if (log_fieldTol) then
    if (fieldTol.ge.0.0d0) then
        write(iow,'(3X,A40,E16.9)')adjl("Tolerance on field error:",40), fieldTol
        write(6  ,'(3X,A40,E16.9)')adjl("Tolerance on field error:",40), fieldTol
    else
        write(ERROR_MESSAGE,'("Tolerance on field error is negative:",E16.9)') fieldTol
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    fieldTol = dflt_fieldTol
    write(iow,'(3X,A40)')adjl("Tolerance on field error not found",40)
    write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), fieldTol
    write(6  ,'(3X,A40)')adjl("Tolerance on field error not found.",40)
    write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), fieldTol
endif


if (log_freeEnergyTol) then
    if (freeEnergyTol.ge.0.0d0) then
        write(iow,'(3X,A40,E16.9)')adjl("Tolerance on free energy error:",40), freeEnergyTol
        write(6  ,'(3X,A40,E16.9)')adjl("Tolerance on free energy error:",40), freeEnergyTol
    else
        write(ERROR_MESSAGE,'("Tolerance on free energy error is negative:",E16.9)') freeEnergyTol
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    freeEnergyTol = dflt_freeEnergyTol
    write(iow,'(3X,A40)')adjl("Tolerance on free energy error not found",40)
    write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), freeEnergyTol
    write(6  ,'(3X,A40)')adjl("Tolerance on free energy error not found.",40)
    write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), freeEnergyTol
endif


if (log_freeEnergyTolForDelta) then
    if (freeEnergyTolForDelta.ge.0.0d0) then
        write(iow,'(3X,A40,E16.9)')adjl("Tolerance on energy error for delta:",40), freeEnergyTolForDelta
        write(6  ,'(3X,A40,E16.9)')adjl("Tolerance on energy error for delta:",40), freeEnergyTolForDelta
    else
        write(ERROR_MESSAGE,'("Tolerance on energy error for delta is negative:",E16.9)') freeEnergyTolForDelta
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    freeEnergyTolForDelta = dflt_freeEnergyTolForDelta
    write(iow,'(3X,A40)')adjl("Tol on energy error for delta not found",40)
    write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), freeEnergyTolForDelta
    write(6  ,'(3X,A40)')adjl("Tol on energy error for delta not found.",40)
    write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40), freeEnergyTolForDelta
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-----------------------------------OUTPUT FREQUENCY-----------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------OUTPUT FREQUENCY-----------------------------------',85)


if (log_exportPhiGeneral) then
    if (exportPhiGeneral.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export density profiles:",40), exportPhiGeneral
        write(6  ,'(3X,A40,I9)')adjl("Export density profiles:",40), exportPhiGeneral
    else
        exportPhiGeneral = 0
        write(iow,'(3X,A40,I9)')adjl("Export density profiles:",40), exportPhiGeneral
        write(6  ,'(3X,A40,I9)')adjl("Export density profiles:",40), exportPhiGeneral
    endif
else
    exportPhiGeneral = dflt_exportPhiGeneral
    write(iow,'(3X,A40,I9)')adjl("Export density profiles:",40), exportPhiGeneral
    write(6  ,'(3X,A40,I9)')adjl("Export density profiles:",40), exportPhiGeneral
endif


if (log_exportPhiIndividual) then
    if (exportPhiIndividual.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export individual density profiles:",40), exportPhiIndividual
        write(6  ,'(3X,A40,I9)')adjl("Export individual density profiles:",40), exportPhiIndividual
    else
        exportPhiIndividual = 0
        write(iow,'(3X,A40,I9)')adjl("Export individual density profiles:",40), exportPhiIndividual
        write(6  ,'(3X,A40,I9)')adjl("Export individual density profiles:",40), exportPhiIndividual
    endif
else
    exportPhiIndividual = dflt_exportPhiIndividual
    write(iow,'(3X,A40,I9)')adjl("Export individual density profiles:",40), exportPhiIndividual
    write(6  ,'(3X,A40,I9)')adjl("Export individual density profiles:",40), exportPhiIndividual
endif


if (log_exportField) then
    if (exportField.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export field:",40), exportField
        write(6  ,'(3X,A40,I9)')adjl("Export field:",40), exportField
    else
        exportField = 0
        write(iow,'(3X,A40,I9)')adjl("Export field:",40), exportField
        write(6  ,'(3X,A40,I9)')adjl("Export field:",40), exportField
    endif
else
    exportField = dflt_exportField
    write(iow,'(3X,A40,I9)')adjl("Export field:",40), exportField
    write(6  ,'(3X,A40,I9)')adjl("Export field:",40), exportField
endif


if (log_exportFieldBinary) then
    if (exportFieldBinary.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export binary field:",40), exportFieldBinary
        write(6  ,'(3X,A40,I9)')adjl("Export binary field:",40), exportFieldBinary
    else
        exportFieldBinary = 0
        write(iow,'(3X,A40,I9)')adjl("Export binary field:",40), exportFieldBinary
        write(6  ,'(3X,A40,I9)')adjl("Export binary field:",40), exportFieldBinary
    endif
else
    exportFieldBinary = dflt_exportFieldBinary
    write(iow,'(3X,A40,I9)')adjl("Export binary field:",40), exportFieldBinary
    write(6  ,'(3X,A40,I9)')adjl("Export binary field:",40), exportFieldBinary
endif


if (log_exportPropagators) then
    if (exportPropagators.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export propagators:",40), exportPropagators
        write(6  ,'(3X,A40,I9)')adjl("Export propagators:",40), exportPropagators
    else
        exportPropagators = 0
        write(iow,'(3X,A40,I9)')adjl("Export propagators:",40), exportPropagators
        write(6  ,'(3X,A40,I9)')adjl("Export propagators:",40), exportPropagators
    endif
else
    exportPropagators = dflt_exportPropagators
    write(iow,'(3X,A40,I9)')adjl("Export propagators:",40), exportPropagators
    write(6  ,'(3X,A40,I9)')adjl("Export propagators:",40), exportPropagators
endif


if (log_exportBrushThickness) then
    if (exportBrushThickness.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export brush thickness:",40), exportBrushThickness
        write(6  ,'(3X,A40,I9)')adjl("Export brush thickness:",40), exportBrushThickness
    else
        exportBrushThickness = 0
        write(iow,'(3X,A40,I9)')adjl("Export brush thickness:",40), exportBrushThickness
        write(6  ,'(3X,A40,I9)')adjl("Export brush thickness:",40), exportBrushThickness
    endif
else
    exportBrushThickness = dflt_exportBrushThickness
    write(iow,'(3X,A40,I9)')adjl("Export brush thickness:",40), exportBrushThickness
    write(6  ,'(3X,A40,I9)')adjl("Export brush thickness:",40), exportBrushThickness
endif


if (log_exportChainsPerArea) then
    if (exportChainsPerArea.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export chains per area profiles:",40), exportChainsPerArea
        write(6  ,'(3X,A40,I9)')adjl("Export chains per area profiles:",40), exportChainsPerArea
    else
        exportChainsPerArea = 0
        write(iow,'(3X,A40,I9)')adjl("Export chains per area profiles:",40), exportChainsPerArea
        write(6  ,'(3X,A40,I9)')adjl("Export chains per area profiles:",40), exportChainsPerArea
    endif
else
    exportChainsPerArea = dflt_exportChainsPerArea
    write(iow,'(3X,A40,I9)')adjl("Export chains per area profiles:",40), exportChainsPerArea
    write(6  ,'(3X,A40,I9)')adjl("Export chains per area profiles:",40), exportChainsPerArea
endif


if (log_exportAdsorbedFree) then
    if (exportAdsorbedFree.ge.1) then
        write(iow,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40), exportAdsorbedFree
        write(6  ,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40), exportAdsorbedFree
    else
        exportAdsorbedFree = 0
        write(iow,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40), exportAdsorbedFree
        write(6  ,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40), exportAdsorbedFree
    endif
else
    exportAdsorbedFree = dflt_exportAdsorbedFree
    write(iow,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40), exportAdsorbedFree
    write(6  ,'(3X,A40,I9)')adjl("Export ads vs free density profiles:",40), exportAdsorbedFree
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('----------------------------------BOUNDARY CONDITIONS---------------------------------',85)
write(*  ,'(A85)')adjl('----------------------------------BOUNDARY CONDITIONS---------------------------------',85)


if (log_numDirichletFaces.and.numDirichletFaces>0) then
    write(iow,'(3X,A40,I9)')adjl("Number of Dirichlet (q=0) faces:",40), numDirichletFaces
    write(iow,'(3X,A39)',advance='no') "Face ids:                               "
    write(6  ,'(3X,A40,I9)')adjl("Number of Dirichlet (q=0) faces:",40), numDirichletFaces
    write(6  ,'(3X,A39)',advance='no') "Face ids:                               "
    do ii = 1, numDirichletFaces
        write(iow,'(I3)',advance='no') dirichletFaceId(ii)
        write(6  ,'(I3)',advance='no') dirichletFaceId(ii)
    enddo
    write(iow  ,'(3X,A39)',advance='yes') "Face values:                            "
    write(6  ,'(3X,A39)',advance='yes')   "Face values:                            "
    do ii = 1, numDirichletFaces
        write(iow,'(F3.1)',advance='no') dirichletFaceValue(ii)
        write(6  ,'(F3.1)',advance='no') dirichletFaceValue(ii)
    enddo
    write(iow,*)
    write(6  ,*)
else
    numDirichletFaces = 0
    allocate(dirichletFaceId(1))
    dirichletFaceId(1)=-1
    write(iow,'(3X,A40)')adjl("There are no Dirichlet (q=0) faces.",40)
    write(6  ,'(3X,A40)')adjl("There are no Dirichlet (q=0) faces.",40)
endif


if (log_numNanoparticleFaces.and.numNanoparticleFaces>0) then
    write(iow,'(3X,A40,I9)')adjl("Number of Nanoparticle (q=0) faces:",40), numNanoparticleFaces
    write(iow,'(3X,A39)',advance='no') "Face ids:                               "
    write(6  ,'(3X,A40,I9)')adjl("Number of Nanoparticle (q=0) faces:",40), numNanoparticleFaces
    write(6  ,'(3X,A39)',advance='no') "Face ids:                               "
    do ii = 1, numNanoparticleFaces
        write(iow,'(I3)',advance='no') nanoparticleFaceId(ii)
        write(6  ,'(I3)',advance='no') nanoparticleFaceId(ii)
    enddo
    write(iow,*)
    write(6  ,*)
else
    numNanoparticleFaces = 0
    allocate(nanoparticleFaceId(1))
    nanoparticleFaceId(1)=-1
    write(iow,'(3X,A38)')adjl("There are no nanoparticle (q=0) faces.",40)
    write(6  ,'(3X,A38)')adjl("There are no nanoparticle (q=0) faces.",40)
endif


if (log_periodicity.and.periodicity>0) then
    domainIsPeriodic = .True.

    if (periodicAxisId(1)) then
        write(iow,'(3X,A40)')adjl("The domain is periodic along the x-axis",40)
        write(iow,'(6X,A10,2I3)') "Face ids: ", periodicFaceId(1), periodicFaceId(2)
        write(6  ,'(3X,A40)')adjl("The domain is periodic along the x-axis",40)
        write(6  ,'(6X,A10,2I3)') "Face ids: ", periodicFaceId(1), periodicFaceId(2)
    endif
    if (periodicAxisId(2)) then
        write(iow,'(3X,A40)')adjl("The domain is periodic along the y-axis",40)
        write(iow,'(6X,A10,2I3)') "Face ids: ", periodicFaceId(3), periodicFaceId(4)
        write(6  ,'(3X,A40)')adjl("The domain is periodic along the y-axis",40)
        write(6  ,'(6X,A10,2I3)') "Face ids: ", periodicFaceId(3), periodicFaceId(4)
    endif
    if (periodicAxisId(3)) then
        write(iow,'(3X,A40)')adjl("The domain is periodic along the z-axis",40)
        write(iow,'(6X,A10,2I3)') "Face ids: ", periodicFaceId(5), periodicFaceId(6)
        write(6  ,'(3X,A40)')adjl("The domain is periodic along the z-axis",40)
        write(6  ,'(6X,A10,2I3)') "Face ids: ", periodicFaceId(5), periodicFaceId(6)
    endif
else
    domainIsPeriodic = dflt_domainIsPeriodic
    write(iow,'(3X,A38)')adjl("The domain is not periodic.",40)
    write(6  ,'(3X,A38)')adjl("The domain is not periodic.",40)
endif



write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('---------------------------------HAMAKER INTERACTIONS---------------------------------',85)
write(*  ,'(A85)')adjl('---------------------------------HAMAKER INTERACTIONS---------------------------------',85)


if (log_sigma_polymer) then
    if (sigma_pol>0.0d0) then
        write(iow,'(3X,A40,E16.9,A11)')adjl("Sigma of polymer:",40), sigma_pol, " [Angstrom]"
        write(6  ,'(3X,A40,E16.9,A11)')adjl("Sigma of polymer:",40), sigma_pol, " [Angstrom]"
    else
        write(ERROR_MESSAGE,'("sigma_polymer is negative: ",E16.9," Angstroms")') sigma_pol
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="sigma_polymer not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif


if (log_Hamaker_constant_of_polymer) then
    write(iow,'(3X,A40,E16.9,A10)')adjl("Hamaker constant of polymer:",40), A_pol, " [10-20 J]"
    write(6  ,'(3X,A40,E16.9,A10)')adjl("Hamaker constant of polymer:",40), A_pol, " [10-20 J]"
    A_pol = A_pol * 1.0d-20
else
    A_pol = dflt_polymer_hamaker_constant
    write(iow,'(3X,A40)')adjl("Hamaker constant for pol not found.",40)
    write(iow,'(3X,A40,E16.9,A10)')adjl("It was set to the default value: ",40), A_pol, " [10-20 J]"
    write(6  ,'(3X,A40)')adjl("Hamaker constant for pol not found.",40)
    write(6  ,'(3X,A40,E16.9,A10)')adjl("It was set to the default value: ",40), A_pol, " [10-20 J]"
endif


if (log_wall_distance) then
    write(iow,'(3X,A40,E16.9,A11)')adjl("Wall distance for Hamaker:",40), wall_distance, " [Angstrom]"
    write(6  ,'(3X,A40,E16.9,A11)')adjl("Wall distance for Hamaker:",40), wall_distance, " [Angstrom]"
else
    wall_distance = dflt_wall_distance
    write(iow,'(3X,A24)')adjl("Wall distance not found.",40)
    write(iow,'(3X,A40,E16.9,A11)')adjl("It was set to the default value: ",40), wall_distance, " [Angstrom]"
    write(6  ,'(3X,A24)')adjl("Wall distance not found.",40)
    write(6  ,'(3X,A40,E16.9,A11)')adjl("It was set to the default value: ",40), wall_distance, " [Angstrom]"
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-------------------------------NONBONDED INTERACTIONS---------------------------------',85)
write(*  ,'(A85)')adjl('-------------------------------NONBONDED INTERACTIONS---------------------------------',85)


if (log_eos_type) then
    if (eos_type.eq.eos_helfand) then
        write(iow,'(3X,A40,I9)')adjl("Equation of state: Helfand",40), eos_type
        write(*  ,'(3X,A40,I9)')adjl("Equation of state: Helfand",40), eos_type
    elseif (eos_type.eq.eos_sl) then
        write(iow,'(3X,A45,I9)')adjl("Equation of state: Sanchez-Lacombe",40), eos_type
        write(*  ,'(3X,A45,I9)')adjl("Equation of state: Sanchez-Lacombe",40), eos_type
    else
        write(*  ,'(A45,I11)') 'EOS flag different than 0 (HF) or 1 (SL)', eos_type
        STOP
    endif
    call init_scf_params()
else
    write(iow,'(A45)') "EOS flag not set"
    write(*  ,'(A45)') "EOS flag not set"
    STOP
endif


if (log_eos_coeffs) then
    if (eos_type.eq.eos_helfand) then
        write(iow,'(3X,A40,E16.9,A8)')adjl("Helfand isothermal compressibility:",40), hlf_kappa_T, " [Pa^-1]"
        write(*  ,'(3X,A40,E16.9,A8)')adjl("Helfand isothermal compressibility:",40), hlf_kappa_T, " [Pa^-1]"
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
    squareGradient = .True.
    write(iow,'(3X,A40,E16.9,A14)')adjl("Influence parameter:",45), k_gr_tilde, " [J*m^5/mol^2]"
    write(*  ,'(3X,A40,E16.9,A14)')adjl("Influence parameter:",45), k_gr_tilde, " [J*m^5/mol^2]"
else
    k_gr_tilde      = 0.0d0
    squareGradient = dflt_squareGradient
    write(iow,'(3X,A40,E16.9,A14)')adjl("Influence parameter not found. Auto:",45), k_gr_tilde, " [J*m^5/mol^2]"
    write(*  ,'(3X,A40,E16.9,A14)')adjl("Influence parameter not found. Auto:",45), k_gr_tilde, " [J*m^5/mol^2]"
endif


write(iow,*)
write(*,*)
write(iow,'(A85)')adjl('-------------------------------------SPATIAL MESH------------------------------------',85)
write(*  ,'(A85)')adjl('-------------------------------------SPATIAL MESH------------------------------------',85)


if (log_meshFile) then
    inquire(file=meshFile, exist=file_exists)
    if (.not.file_exists) then
        write(ERROR_MESSAGE,'("Mesh file ",A16," does not exist!")') meshFile
        call exit_with_error(1,1,1,ERROR_MESSAGE)
        STOP
    endif
    write(iow,'(3X,A40,A16)')adjl("Reading mesh from file:",40), meshFile
    write(6  ,'(3X,A40,A16)')adjl("Reading mesh from file:",40), meshFile
else
    meshFile = dflt_meshFile
    write(iow,'(3X,A40)')adjl("Mesh input file not specified.",40)
    write(iow,'(3X,A40,A16)')adjl("Reading default mesh file:",40), meshFile
    write(6  ,'(3X,A40)')adjl("Mesh input file not specified.",40)
    write(6  ,'(3X,A40,A16)')adjl("Reading default mesh file:",40), meshFile

    inquire(file=meshFile, exist=file_exists)
    if (.not.file_exists) then
        write(ERROR_MESSAGE,'("Default mesh file ",A16," does not exist!")') meshFile
        call exit_with_error(1,1,1,ERROR_MESSAGE)
        STOP
    endif
endif

return
!--------------------------------------------------------------------------------!
end subroutine parser_input
