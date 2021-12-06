subroutine parser
!--------------------------------------------------------------------------------!
use parser_vars
use write_helper
use error_handing
use iofiles
!--------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------!
character(100) :: line

integer :: Reason
integer :: i1, id

logical :: log_set_initial_iteration       = .true.
logical :: log_number_of_iterations        = .false.
logical :: log_field_init_scheme           = .false.
logical :: log_output_every                = .false.
logical :: log_profile_dimension           = .false.
logical :: log_use_grafted                 = .false.
logical :: log_n_matrix_seg                = .false.
logical :: log_n_gr_seg                    = .false.
logical :: log_mumps_matrix_type           = .false.
logical :: log_time_integration_scheme     = .false.
logical :: log_grafted_ic_from_delta       = .false.
logical :: log_n_dirichlet_faces           = .false.
logical :: log_n_nanopart_faces            = .false.

logical :: log_temperature                 = .false.
logical :: log_isothermal_compressibility  = .false.
logical :: log_mass_density                = .false.
logical :: log_monomer_mass                = .false.
logical :: log_Rg2_per_mon_matrix          = .false.
logical :: log_Rg2_per_mon_gr              = .false.
logical :: log_chainlen_matrix             = .false.
logical :: log_chainlen_gr                 = .false.
logical :: log_interf_area                 = .false.
logical :: log_sigma_polymer               = .false.
logical :: log_Hamaker_constant_of_polymer = .false.
logical :: log_wall_distance               = .false.
logical :: log_fraction_of_new_field       = .false.
logical :: log_maximum_error               = .false.
!--------------------------------------------------------------------------------!
!parse input file to retrieve simulation parameters
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
        if (index(line,"# temperature") > 0) then
            read(line,*) temp
            log_temperature = .true.
        elseif (index(line,"# isothermal compressibility") > 0) then
            read(line,*) kappa_T
            log_isothermal_compressibility = .true.
        elseif (index(line,"# mass density") > 0) then
            read(line,*) massden
            log_mass_density = .true.
        elseif (index(line,"# monomer mass") > 0) then
            read(line,*) mon_mass
            log_monomer_mass = .true.
        elseif (index(line,"# Rg2 per monomer for matrix chains") > 0) then
            read(line,*) Rg2_per_mon_matrix
            log_Rg2_per_mon_matrix = .true.
        elseif (index(line,"# Rg2 per monomer for grafted chains") > 0) then
            read(line,*) Rg2_per_mon_gr
            log_Rg2_per_mon_gr = .true.
        elseif (index(line,"# chain length of matrix chains") > 0) then
            read(line,*) chainlen_matrix
            log_chainlen_matrix = .true.
        elseif (index(line,"# chain length of grafted chains ") > 0) then
            read(line,*) chainlen_gr
            log_chainlen_gr = .true.
        elseif (index(line,"# interfacial area") > 0) then
            read(line,*) interf_area
            log_interf_area = .true.
        elseif (index(line,"# sigma polymer") > 0) then
            read(line,*) sigma_pol
            log_sigma_polymer = .true.
        elseif (index(line,"# Hamaker constant of polymer") > 0) then
            read(line,*) A_pol
            log_Hamaker_constant_of_polymer = .true.
        elseif (index(line,"# wall distance") > 0) then
            read(line,*) wall_distance
            log_wall_distance = .true.
        elseif (index(line,"# fraction of new field") > 0) then
            read(line,*) frac
            log_fraction_of_new_field = .true.
        elseif (index(line,"# maximum error") > 0) then
            read(line,*) max_error_tol
            log_maximum_error = .true.
        elseif (index(line,"# set initial iteration") > 0) then
            read(line,*) init_iter
            log_set_initial_iteration= .true.
        elseif (index(line,"# number of iterations") > 0) then
            read(line,*) iterations
            log_number_of_iterations = .true.
        elseif (index(line,"# initialize field") > 0) then
            read(line,*) field_init_scheme
            log_field_init_scheme= .true.
        elseif (index(line,"# output every") > 0) then
            read(line,*) output_every
            log_output_every = .true.
        elseif (index(line,"# profile dimension") > 0) then
            read(line,*) prof_dim
            log_profile_dimension = .true.
        elseif (index(line,"# use grafted") > 0) then
            read(line,*) use_grafted
            log_use_grafted = .true.
        elseif (index(line,"# number of matrix segments") > 0) then
            read(line,*) ns_matrix_ed, ns_matrix_conv
            log_n_matrix_seg = .true.
        elseif (index(line,"# number of grafted segments") > 0) then
            read(line,*) ns_gr_ed, ns_gr_conv
            log_n_gr_seg = .true.
        elseif (index(line,"# mumps matrix type") > 0) then
            read(line,*) mumps_matrix_type
            log_mumps_matrix_type = .true.
        elseif (index(line,"# time integration scheme") > 0) then
            read(line,*) time_integration_scheme
            log_time_integration_scheme = .true.
        elseif (index(line,"# calculate grafted initial condition using delta function") > 0) then
            read(line,*) grafted_ic_from_delta
            log_grafted_ic_from_delta = .true.
        elseif (index(line,"# n dirichlet faces") > 0) then
            read(line,*) n_dirichlet_faces
            if (n_dirichlet_faces > 0) then
                allocate(ids_dirichlet_faces(n_dirichlet_faces))
                allocate(sigma_plate(n_dirichlet_faces))
                allocate(A_plate(n_dirichlet_faces))
                do i1 = 1, n_dirichlet_faces
                    read(256,*) id, sigma_plate(i1), A_plate(i1)
                    ids_dirichlet_faces(i1) = id
                enddo
                A_plate = A_plate * 1e-20
                log_n_dirichlet_faces = .true.
            endif
        elseif (index(line,"# n nanoparticle faces") > 0) then
            read(line,*) n_nanopart_faces
            if (n_nanopart_faces > 0) then
                allocate(ids_nanopart_faces(n_nanopart_faces))
                allocate(center_np(3,n_nanopart_faces))
                allocate(radius_np(n_nanopart_faces))
                allocate(sigma_np(n_nanopart_faces))
                allocate(A_np(n_nanopart_faces))
                do i1 = 1, n_nanopart_faces
                    read(256,*) id, radius_np(i1), center_np(1,i1), center_np(2,i1), center_np(3,i1), &
                                                          & sigma_np(i1), A_np(i1)
                    ids_nanopart_faces(i1) = id
                enddo
                A_np = A_np * 1e-20
                log_n_nanopart_faces = .true.
            endif
        endif
   endif
enddo

close(256)

!check parameters of input file
if (log_set_initial_iteration) then
    if (init_iter.lt.0) then
        write(ERROR_MESSAGE,'("Wrong value of initial iteration.",I16)') init_iter
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    init_iter = 0
endif
if (init_iter.eq.0) then
    write(iow,'(3X,A40,I16)')adjl("Fresh simulation starting from iter:",40), init_iter
    write(6  ,'(3X,A40,I16)')adjl("Fresh simulation starting from iter:",40), init_iter
elseif (init_iter.gt.0) then
    write(iow,'(3X,A40,I16)')adjl("Simulation restarting from iter:",40), init_iter
    write(6  ,'(3X,A40,I16)')adjl("Simulation restarting from iter:",40), init_iter
endif

if (log_number_of_iterations) then
    if (iterations>0) then
        write(iow,'(3X,A40,I16)')adjl("Maximum number of iterations:",40),iterations
        write(6  ,'(3X,A40,I16)')adjl("Maximum number of iterations:",40),iterations
    else
        write(ERROR_MESSAGE,'("Maximum number of iterations is negative:",I10)') iterations
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    iterations = 1
    write(iow,'(3X,A40)')adjl("Max number of iter not found.",40)
    write(iow,'(3X,A40,I16)')adjl("It was set to the default value:",40),iterations
    write(6  ,'(3X,A40)')adjl("Max number of iter not found.",40)
    write(6  ,'(3X,A40,I16)')adjl("It was set to the default value:",40),iterations
endif

if (log_field_init_scheme) then
    if (field_init_scheme==0) then
        write(iow,'(3X,A34)')adjl("Field will be initialized to zero.",40)
        write(6  ,'(3X,A34)')adjl("Field will be initialized to zero.",40)
    elseif (field_init_scheme==1) then
        write(iow,'(A43,A15)')adjl("Field will be read from file:",40),field_in_filename
        write(6  ,'(A43,A15)')adjl("Field will be read from file:",40),field_in_filename
    elseif (field_init_scheme==2) then
        write(iow,'(A43)')adjl("Field: -kapa at Dir. BCs and 0 elsewhere.",43)
        write(6  ,'(A43)')adjl("Field: -kapa at Dir. BCs and 0 elsewhere.",43)
    else
        write(ERROR_MESSAGE,'("Incorrect field initialization value. Choose between 1-3.",I16)') init_iter
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
endif
if (.not.log_field_init_scheme) then
    write(iow,'(A40)')adjl("Field will be initialized to zero.",40)
    write(6  ,'(A40)')adjl("Field will be initialized to zero.",40)
endif

if (log_output_every) then
    if (output_every.ge.1) then
        write(iow,'(3X,A40,1x,I15)')adjl("Data is dumped every (steps):",40),output_every
        write(6  ,'(3X,A40,1x,I15)')adjl("Data is dumped every (steps):",40),output_every
    else
        output_every = 1
        write(iow,'(3X,A40,1X,I15)')adjl("Data is dumped every (steps):",40),output_every
        write(6  ,'(3X,A40,1X,I15)')adjl("Data is dumped every (steps):",40),output_every
    endif
else
    output_every = 1
endif

if (log_profile_dimension) then
    if ( (prof_dim.ge.1).and.(prof_dim.le.3) ) then
        write(iow,'(3X,A40,I16)')adjl("Profile dimension:",40),prof_dim
        write(6  ,'(3X,A40,I16)')adjl("Profile dimension:",40),prof_dim
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

if (log_use_grafted.and.use_grafted.ge.1) then
    use_grafted = 1
    write(iow,'(3X,A40,1x,I15)')adjl("The system includes grafted chains:",40),use_grafted
    write(6  ,'(3X,A40,1x,I15)')adjl("The system includes grafted chains:",40),use_grafted
else
    use_grafted = 0
    write(iow,'(3X,A40,1x,I15)')adjl("System does not include grafted chains:",40),use_grafted
    write(6  ,'(3X,A40,1x,I15)')adjl("System does not include grafted chains:",40),use_grafted
endif

if (log_n_matrix_seg) then
    if (ns_matrix_ed>0 .and. ns_matrix_conv>0) then
        write(iow,'(3X,A40,I9,I7)')adjl("Number of matrix segments:",40),ns_matrix_ed, ns_matrix_conv
        write(6  ,'(3X,A40,I9,I7)')adjl("Number of matrix segments:",40),ns_matrix_ed, ns_matrix_conv
        if (mod(ns_matrix_ed,2).ne.0 .or. mod(ns_matrix_conv,2).ne.0) then
            write(ERROR_MESSAGE,'("ns_matrix is not an even number: ",I16,I16)') ns_matrix_ed, ns_matrix_conv
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        write(ERROR_MESSAGE,'("ns_matrix is negative: ",I16,I16)') ns_matrix_ed, ns_matrix_conv
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Number of matrix segments was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (use_grafted.eq.1) then
    if (log_n_gr_seg) then
        if (ns_gr_ed>0 .and. ns_gr_conv>0) then
            write(iow,'(3X,A40,I9,I7)')adjl("Number of grafted segments:",40),ns_gr_ed, ns_gr_conv
            write(6  ,'(3X,A40,I9,I7)')adjl("Number of grafted segments:",40),ns_gr_ed, ns_gr_conv
            if (mod(ns_gr_ed,2).ne.0 .or. mod(ns_gr_conv,2).ne.0) then
                write(ERROR_MESSAGE,'("ns_grafted is not an even number: ",I16,I16)') ns_gr_ed, ns_gr_conv
                call exit_with_error(1,1,1,ERROR_MESSAGE)
            endif
        else
            write(ERROR_MESSAGE,'("ns_grafted is negative: ",I16,I16)') ns_gr_ed, ns_gr_conv
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Number of grafted segments was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_grafted_ic_from_delta) then
        write(iow,'(3X,A40,1X,I15)')adjl("Grafted ic from delta:",40),grafted_ic_from_delta
        write(6  ,'(3X,A40,1X,I15)')adjl("Grafted ic from delta:",40),grafted_ic_from_delta
    else
        ERROR_MESSAGE="Grafted ic from delta variable was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
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

if (log_time_integration_scheme) then
    if (time_integration_scheme.eq.1) then
        write(iow,'(3X,A40,I16)')adjl("Time integr with uniform spacing:",40),time_integration_scheme
        write(6  ,'(3X,A40,I16)')adjl("Time integr with uniform spacing:",40),time_integration_scheme
    elseif ( time_integration_scheme.eq.2) then
        write(iow,'(3X,A40,I16)')adjl("Time integr with non-uniform spacing:",40),time_integration_scheme
        write(6  ,'(3X,A40,I16)')adjl("Time integr with non-uniform spacing:",40),time_integration_scheme
    else
        write(ERROR_MESSAGE,'(''Time integration scheme does not exist:'',I16)') time_integration_scheme
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    time_integration_scheme = 1
    write(iow,'(3X,A40)')adjl("Time integration scheme not found.",40)
    write(iow,'(3X,A40,I8)')adjl("Integration with Simpson Rule:",40),time_integration_scheme
    write(6  ,'(3X,A40)')adjl("Time integration scheme not found.",40)
    write(6  ,'(3X,A40,I8)')adjl("Integration with Simpson Rule:",40),time_integration_scheme
endif

if (log_n_dirichlet_faces.and.n_dirichlet_faces>0) then
    write(iow,'(3X,A40,I16)')adjl("Number of Dirichlet (q=0) faces:",40), n_dirichlet_faces
    write(iow,'(3X,A39)',advance='no')"Face ids:                               "
    write(6  ,'(3X,A40,I16)')adjl("Number of Dirichlet (q=0) faces:",40), n_dirichlet_faces
    write(6  ,'(3X,A39)',advance='no')"Face ids:                               "
    do i1 = 1, n_dirichlet_faces
        write(iow,'(I3)',advance='no') ids_dirichlet_faces(i1)
        write(6  ,'(I3)',advance='no') ids_dirichlet_faces(i1)
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
    write(iow,'(3X,A40,I8)')adjl("Number of Nanoparticle (q=0) faces:",40), n_nanopart_faces
    write(iow,'(3X,A39)',advance='no')"Face ids:                               "
    write(6  ,'(3X,A40,I8)')adjl("Number of Nanoparticle (q=0) faces:",40), n_nanopart_faces
    write(6  ,'(3X,A39)',advance='no')"Face ids:                               "
    do i1 = 1, n_nanopart_faces
        write(iow,'(I3)',advance='no') ids_nanopart_faces(i1)
        write(6  ,'(I3)',advance='no') ids_nanopart_faces(i1)
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

if (log_temperature) then
    if (temp>0) then
        write(iow,'(3X,A40,E16.9,A4)')adjl("Temperature:",40),temp," [K]"
        write(6  ,'(3X,A40,E16.9,A4)')adjl("Temperature:",40),temp," [K]"
    else
        write(ERROR_MESSAGE,'("Temperature is negative: ",E16.9," K")') temp
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Temperature was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_isothermal_compressibility) then
    if (kappa_T>0) then
        write(iow,'(3X,A40,E16.9,A12)')adjl("Isothermal compressibility:",40),kappa_T*1.e9," [10-9 Pa-1]"
        write(6  ,'(3X,A40,E16.9,A12)')adjl("Isothermal compressibility:",40),kappa_T*1.e9," [10-9 Pa-1]"
    else
        write(ERROR_MESSAGE,'("Isothermal compressibility is negative: ",E16.9," x10-9 Pa-1")') kappa_T*1.e9
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Isothermal compressibility was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

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
        write(iow,'(3X,A40,E16.9,A9)')adjl("Monomer mass:",40),mon_mass,"[g/mol]"
        write(6  ,'(3X,A40,E16.9,A9)')adjl("Monomer mass:",40),mon_mass,"[g/mol]"
    else
        write(ERROR_MESSAGE,'("Monomer mass is negative: ",E16.9)') mon_mass
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Monomer mass not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_Rg2_per_mon_matrix) then
    if (Rg2_per_mon_matrix>0) then
        write(iow,'(3X,A40,E16.9)')adjl("Rg2 per matrix monomer:",40),Rg2_per_mon_matrix
        write(6  ,'(3X,A40,E16.9)')adjl("Rg2 per matrix monomer:",40),Rg2_per_mon_matrix
    else
        write(ERROR_MESSAGE,'("Rg2 per matrix monomer is negative: ",E16.9)') Rg2_per_mon_matrix
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Rg2 per matrix monomer was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_chainlen_matrix) then
    if (chainlen_matrix>0) then
        write(iow,'(3X,A40,E16.9)')adjl("Chain length of matrix chains:",40),chainlen_matrix
        write(6  ,'(3X,A40,E16.9)')adjl("Chain length of matrix chains:",40),chainlen_matrix
    else
        write(ERROR_MESSAGE,'("Chain length of matrix chains is negative: ",E16.9)') chainlen_matrix
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Chain length of matrix chains was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (use_grafted.eq.1) then
    if (log_Rg2_per_mon_gr) then
        if ( Rg2_per_mon_gr > 0 ) then
            write(iow,'(3X,A40,E16.9)')adjl("Rg2 per grafted monomer:",40),Rg2_per_mon_gr
            write(6  ,'(3X,A40,E16.9)')adjl("Rg2 per grafted monomer:",40),Rg2_per_mon_gr
        else
            write(ERROR_MESSAGE,'("Rg2 per grafted monomer is negative: ",E16.9)') Rg2_per_mon_gr
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Rg2 per matrix monomer was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    if (log_chainlen_gr) then
        if ( chainlen_gr > 0 ) then
            write(iow,'(3X,A40,E16.9)')adjl("Chain length of grafted chains:",40),chainlen_gr
            write(6  ,'(3X,A40,E16.9)')adjl("Chain length of grafted chains:",40),chainlen_gr
        else
            write(ERROR_MESSAGE,'("Chain length of grafted chains is negative: ",E16.9)') chainlen_gr
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        ERROR_MESSAGE="Chain length of grafted chains was not detected."
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
endif

if (log_interf_area) then
    if (interf_area>0) then
        write(iow,'(3X,A40,E16.9,A13)')adjl("Interface area:",40),interf_area," [Angstrom^2]"
        write(6  ,'(3X,A40,E16.9,A13)')adjl("Interface area:",40),interf_area," [Angstrom^2]"
    else
        write(ERROR_MESSAGE,'("Interface area is negative: ",E16.9," Angstrom^2")') interf_area
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Interface area was not detected."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

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
        write(iow,'(3X,A40,E16.9)')adjl("Maximum tolerance error:",40),max_error_tol
        write(6  ,'(3X,A40,E16.9)')adjl("Maximum tolerance error:",40),max_error_tol
    else
        write(ERROR_MESSAGE,'("Maximum tolerance error is negative:",E16.9)') max_error_tol
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    max_error_tol = 0.d0
    write(iow,'(3X,A40)')adjl("Max error not found",40)
    write(iow,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),max_error_tol
    write(6  ,'(3X,A40)')adjl("Max error not found.",40)
    write(6  ,'(3X,A40,E16.9)')adjl("It was set to the default value:",40),max_error_tol
endif

return
!--------------------------------------------------------------------------------!
end subroutine parser
