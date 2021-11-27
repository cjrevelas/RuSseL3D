subroutine scfinout
!--------------------------------------------------------------------------------!
use xdata
use constants
use write_helper
use error_handing
!--------------------------------------------------------------------------------!
implicit none

CHARACTER(100) :: line
CHARACTER(15) :: input_filename = 'gaussdat.in.txt'
character(12) :: field_filename = 'field.in.bin'
integer :: Reason

logical :: log_domain_geometry = .false.
logical :: log_timestep = .false.
logical :: log_chain_length = .false.
logical :: log_temperature = .false.
logical :: log_Hamaker_constant_of_polymer = .false.
logical :: log_Hamaker_constant_of_solid = .false.
logical :: log_mass_density = .false.
logical :: log_isothermal_compressibility = .false.
logical :: log_characteristic_ratio = .false.
logical :: log_number_of_iterations = .false.
logical :: log_maximum_error = .false.
logical :: log_fraction_of_new_field = .false.
logical :: log_monomer_mass = .false.
logical :: log_sphere_radius = .false.
logical :: log_sigma_polymer = .false.
logical :: log_sigma_solid = .false.
logical :: log_read_field = .false.
logical :: log_mix_coef_fraction = .false.
logical :: log_mix_coef_kapa = .false.
logical :: log_n_dirichlet_faces = .false.
logical :: log_convergence_scheme = .false.
logical :: log_mumps_matrix_type = .false.

!*******************************************************************!
!                    Read the input parameter file                  !
!*******************************************************************!

write(iow,'(/''*Reading parameters from file'',A16,/)') input_filename
write(6  ,'(/''*Reading parameters from file'',A16,/)') input_filename

INQUIRE(FILE=input_filename, EXIST=FILE_EXISTS)

if (FILE_EXISTS) then
    open(unit=256, file = input_filename)
else
    write(ERROR_MESSAGE,'(''File '',A15,'' does not exist!'')')input_filename
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

do
    read(256,'(A100)',IOSTAT=Reason) line

    if (Reason > 0)  then
        write(*,*)"Something went wrong!"
    elseif (Reason < 0) then
        write(*,*)"Input parameter file was read!"
        exit
    else
        if (index(line,'# domain geometry') > 0) then
            read(line,'(I9)') iseed
            log_domain_geometry = .true.
        elseif (index(line,'# timestep') > 0) then
            read(line,'(I6)') ns
            log_timestep = .true.
        elseif (index(line,'# chain length') > 0) then
            read(line,'(E16.9)') chainlen
            log_chain_length = .true.
        elseif (index(line,'# temperature') > 0) then
            read(line,'(E16.9)') temp
            log_temperature = .true.
        elseif (index(line,'# Hamaker constant of polymer') > 0) then
            read(line,'(E16.9)') Aps
            log_Hamaker_constant_of_polymer = .true.
        elseif (index(line,'# Hamaker constant of solid') > 0) then
            read(line,'(E16.9)') Asio2
            log_Hamaker_constant_of_solid = .true.
        elseif (index(line,'# mass density') > 0) then
            read(line,'(E16.9)') massden
            log_mass_density = .true.
        elseif (index(line,'# isothermal compressibility') > 0) then
            read(line,'(E16.9)') kappa_T
            log_isothermal_compressibility = .true.
        elseif (index(line,'# characteristic ratio') > 0) then
            read(line,'(E16.9)') CN
            log_characteristic_ratio = .true.
        elseif (index(line,'# number of iterations') > 0) then
            read(line,'(I10)') iterations
            log_number_of_iterations = .true.
        elseif (index(line,'# maximum error') > 0) then
            read(line,'(F16.4)') max_error_tol
            log_maximum_error = .true.
        elseif (index(line,'# fraction of new field') > 0) then
            read(line,'(F16.9)') fraction
            log_fraction_of_new_field = .true.
        elseif (index(line,'# monomer mass') > 0) then
            read(line,'(F16.9)') mon_mass
            log_monomer_mass = .true.
        elseif (index(line,'# sphere radius') > 0) then
            read(line,'(F16.9)') sphere_radius
            log_sphere_radius = .true.
        elseif (index(line,'# sigma polymer') > 0) then
            read(line,'(F16.9)') sigma1
            log_sigma_polymer = .true.
        elseif (index(line,'sigma solid') > 0) then
            read(line,'(F16.4)') sigma2
            log_sigma_solid = .true.
        elseif (index(line,'# read field') > 0) then
            read(line,'(I10)') readfield
            log_read_field= .true.
        elseif (index(line,'# mix coef fraction') > 0) then
            read(line,'(F16.9)') mix_coef_frac
            log_mix_coef_fraction = .true.
        elseif (index(line,'# mix coef kapa') > 0) then
            read(line,'(F16.9)') mix_coef_kapa
            log_mix_coef_kapa = .true.
        elseif (index(line,'# convergence scheme') > 0) then
            read(line,'(I10)') scheme_type
            log_convergence_scheme = .true.
        elseif (index(line,'# mumps matrix type') > 0) then
            read(line,'(I10)') mumps_matrix_type
            log_mumps_matrix_type = .true.
        elseif (index(line,'# n dirichlet faces') > 0) then
            read(line,'(I10)') n_dirichlet_faces
            allocate(ids_dirichlet_faces(n_dirichlet_faces))
            do i1 = 1, n_dirichlet_faces
                 read(256,'(I10)')id
                 ids_dirichlet_faces(i1) = id
            enddo
            log_n_dirichlet_faces = .true.
        endif
   endif
enddo

close(256)

!*******************************************************************!
!              Check the inputs of the  parameter file              !
!*******************************************************************!

if (log_domain_geometry) then
    if (iseed.eq.0) then
        write(iow,'(3x,A40,1x,A4)')adjl('Domain geometry:',40),'FILM'
        write(6  ,'(3x,A40,1x,A4)')adjl('Domain geometry:',40),'FILM'
    elseif (iseed.eq.1) then
        write(iow,'(3x,A40,1x,A6)')adjl('Domain geometry:',40),'SPHERE'
        write(6  ,'(3x,A40,1x,A6)')adjl('Domain geometry:',40),'SPHERE'
    else
        ERROR_MESSAGE="Wrong domain geometry! Please choose a value between 0 (film) and 1 (sphere)"
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE="Domain geometry was not detected.."
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_timestep) then
    if ( ns > 0) then
        write(iow,'(3x,A40,I16)')adjl('Contour length discretization:',40),ns
        write(6  ,'(3x,A40,I16)')adjl('Contour length discretization:',40),ns
        if (mod(ns,2).ne.0) then
            write(ERROR_MESSAGE,'(''ns is not an even number: '',I16)') ns
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    else
        write(ERROR_MESSAGE,'(''ns is negative: '',I16)') ns
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='timestep was not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_chain_length) then
    if ( chainlen > 0 ) then
        write(iow,'(3x,A40,E16.9)')adjl('Chain length:',40),chainlen
        write(6  ,'(3x,A40,E16.9)')adjl('Chain length:',40),chainlen
    else
        write(ERROR_MESSAGE,'(''chain length is negative: '',E16.9)') chainlen
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='chain length was not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_temperature) then
    if ( temp > 0) then
        write(iow,'(3x,A40,E16.9,A2)')adjl('temperature:',40),temp,' K'
        write(6  ,'(3x,A40,E16.9,A2)')adjl('temperature:',40),temp,' K'
    else
        write(ERROR_MESSAGE,'(''temperature is negative: '',E16.9,'' K'')') temp
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='temperature was not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_Hamaker_constant_of_polymer) then
    write(iow,'(3x,A40,E16.9,A13)')adjl('Hamaker constant of polymer:',40),Aps,' 10-20 joules'
    write(6  ,'(3x,A40,E16.9,A13)')adjl('Hamaker constant of polymer:',40),Aps,' 10-20 joules'
else
    Aps = 0.d0
    write(iow,'(3x,A40)')adjl('Hamaker constant for pol not found..',40)
    write(iow,'(3x,A40,E16.9,A13)')adjl('---It was set to the default value: ',40),Aps,' 10-20 joules'
    write(6  ,'(3x,A40)')adjl('Hamaker constant for pol not found..',40)
    write(6  ,'(3x,A40,E16.9,A13)')adjl('---It was set to the default value: ',40),Aps,' 10-20 joules'
endif

if (log_Hamaker_constant_of_solid) then
    write(iow,'(3x,A40,E16.9,A13)')adjl('Hamaker constant of solid:',40),Asio2,' 10-20 joules'
    write(6  ,'(3x,A40,E16.9,A13)')adjl('Hamaker constant of solid:',40),Asio2,' 10-20 joules'
else
    Asio2 = 0.d0
    write(iow,'(3x,A40)')adjl('Hamaker constant for solid not found..',40)
    write(iow,'(3x,A40,E16.9,A13)')adjl('---It was set to the default value: ',40),Asio2,' 10-20 joules'
    write(6  ,'(3x,A40)')adjl('Hamaker constant for solid not found..',40)
    write(6  ,'(3x,A40,E16.9,A13)')adjl('---It was set to the default value: ',40),Asio2,' 10-20 joules'
endif

if (log_mass_density) then
    if ( massden > 0) then
        write(iow,'(3x,A40,E16.9,A6)')adjl('Mass density:',40),massden,' g/cm3'
        write(6  ,'(3x,A40,E16.9,A6)')adjl('Mass density:',40),massden,' g/cm3'
    else
        write(ERROR_MESSAGE,'(''Mass density is negative: '',E16.9,'' g/cm3'')') massden
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='Mass density was not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_isothermal_compressibility) then
    if ( kappa_T > 0) then
        write(iow,'(3x,A40,E16.9,A11)')adjl('Isothermal compressibility:',40),kappa_T*1.e9,' x10-9 Pa-1'
        write(6  ,'(3x,A40,E16.9,A11)')adjl('Isothermal compressibility:',40),kappa_T*1.e9,' x10-9 Pa-1'
    else
        write(ERROR_MESSAGE,'(''Isothermal compressibility is negative: '',E16.9,'' x10-9 Pa-1'')') kappa_T*1.e9
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='Isothermal compressibility was not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_characteristic_ratio) then
    if ( CN > 0) then
        write(iow,'(3x,A40,E16.9)')adjl('characteristic_ratio:',40),CN
        write(6  ,'(3x,A40,E16.9)')adjl('characteristic_ratio:',40),CN
    else
        write(ERROR_MESSAGE,'(''characteristic_ratio is negative:'',E16.9)') CN
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='characteristic_ratio was not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_number_of_iterations) then
    if ( iterations > 0) then
        write(iow,'(3x,A40,I16)')adjl('maximum number of iterations:',40),iterations
        write(6  ,'(3x,A40,I16)')adjl('maximum number of iterations:',40),iterations
    else
        write(ERROR_MESSAGE,'(''maximum number of iterations is negative:'',I10)') iterations
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    iterations = 1
    write(iow,'(3x,A40)')adjl('max number of iter not found...',40)
    write(iow,'(3x,A40,I16)')adjl('---It was set to the default value:',40),iterations
    write(6  ,'(3x,A40)')adjl('max number of iter not found...',40)
    write(6  ,'(3x,A40,I16)')adjl('---It was set to the default value:',40),iterations
endif

if (log_maximum_error) then
    if ( max_error_tol.ge.0.d0) then
        write(iow,'(3x,A40,E16.9)')adjl('maximum tolerance error:',40),max_error_tol
        write(6  ,'(3x,A40,E16.9)')adjl('maximum tolerance error:',40),max_error_tol
    else
        write(ERROR_MESSAGE,'(''maximum tolerance error is negative:'',E16.9)') max_error_tol
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    max_error_tol = 0.d0
    write(iow,'(3x,A40)')adjl('max error not found..',40)
    write(iow,'(3x,A40,E16.9)')adjl('---It was set to the default value:',40),max_error_tol
    write(6  ,'(3x,A40)')adjl('max error not found..',40)
    write(6  ,'(3x,A40,E16.9)')adjl('---It was set to the default value:',40),max_error_tol
endif

if (log_fraction_of_new_field) then
    if ( fraction.ge.0 .and. fraction.le.1) then
        write(iow,'(3x,A40,E16.9)')adjl('Initial fraction of new field:',40),fraction
        write(6  ,'(3x,A40,E16.9)')adjl('Initial fraction of new field:',40),fraction
    else
        write(ERROR_MESSAGE,'(''Initial fraction of new field is negative or larger than unity:'',E16.9)') fraction
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    fraction = 1.d0
    write(iow,'(3x,A40)')adjl('No initial fraction of new field..',40)
    write(iow,'(3x,A40,E16.9)')adjl('---It was set to the default value:',40),fraction
    write(6  ,'(3x,A40)')adjl('No initial fraction of new field..',40)
    write(6  ,'(3x,A40,E16.9)')adjl('---It was set to the default value:',40),fraction
endif

if (log_monomer_mass) then
    if ( mon_mass > 0) then
        write(iow,'(3x,A40,E16.9,A7)')adjl('monomer mass:',40),mon_mass,' gr/mol'
        write(6  ,'(3x,A40,E16.9,A7)')adjl('monomer mass:',40),mon_mass,' gr/mol'
    else
        write(ERROR_MESSAGE,'(''monomer mass is negative: '',E16.9)') mon_mass
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='monomer mass not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_sphere_radius) then
    if ( sphere_radius > 0) then
        write(iow,'(3x,A40,E16.9,A10)')adjl('sphere radius:',40),sphere_radius,' Angstroms'
        write(6  ,'(3x,A40,E16.9,A10)')adjl('sphere radius:',40),sphere_radius,' Angstroms'
    else
        write(ERROR_MESSAGE,'(''sphere_radius is negative: '',E16.9)') sphere_radius
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='sphere_radius not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_sigma_polymer) then
    if ( sigma1 > 0) then
        write(iow,'(3x,A40,E16.9,A10)')adjl('sigma of polymer:',40),sigma1,' Angstroms'
        write(6  ,'(3x,A40,E16.9,A10)')adjl('sigma of polymer:',40),sigma1,' Angstroms'
    else
        write(ERROR_MESSAGE,'(''sigma_polymer is negative: '',E16.9,'' Angstroms'')') sigma1
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='sigma_polymer not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_sigma_solid) then
    if ( sigma2 > 0) then
        write(iow,'(3x,A40,E16.9,A10)')adjl('sigma of solid:',40),sigma2,' Angstroms'
        write(6  ,'(3x,A40,E16.9,A10)')adjl('sigma of solid:',40),sigma2,' Angstroms'
    else
        write(ERROR_MESSAGE,'(''sigma_solid is negative: '',E16.9,'' Angstroms'')') sigma2
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    ERROR_MESSAGE='sigma_solid not detected..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

if (log_mix_coef_fraction.and.mix_coef_frac.ge.0.d0) then
    write(iow,'(3x,A40,E16.9)')adjl('mix coef fraction:',40),mix_coef_frac
    write(6  ,'(3x,A40,E16.9)')adjl('mix coef fraction:',40),mix_coef_frac
else
    mix_coef_frac = 1.d0
    write(iow,'(3x,A40)')adjl('mix coef fraction not found or negative..',40)
    write(iow,'(3x,A40,E16.9)')adjl('---It was set to the default value:',40), mix_coef_frac
    write(6  ,'(3x,A40)')adjl('mix coef fraction not found or negative..',40)
    write(6  ,'(3x,A40,E16.9)')adjl('---It was set to the default value:',40), mix_coef_frac
endif

if (log_mix_coef_kapa.and.mix_coef_frac.ge.0.d0) then
    write(iow,'(3x,A40,E16.9)')adjl('mix coef kapa:',40),mix_coef_kapa
    write(6  ,'(3x,A40,E16.9)')adjl('mix coef kapa:',40),mix_coef_kapa
else
    mix_coef_kapa = 10.d0
    write(iow,'(3x,A40)')adjl('mix coef kapa not found or negative..',40)
    write(iow,'(3x,A40,E16.9)')adjl('---It was set to the default value:',40), mix_coef_kapa
    write(6  ,'(3x,A40)')adjl('mix coef kapa not found or negative..',40)
    write(6  ,'(3x,A40,E16.9)')adjl('---It was set to the default value:',40), mix_coef_kapa
endif

if (log_n_dirichlet_faces.and.n_dirichlet_faces>0) then
    write(iow,'(3x,A40,I16)')adjl('Number of Dirichlet (q=0) faces:',40), n_dirichlet_faces
    write(iow,'(3x,A39)',advance='no')'---face ids:                               '
    write(6  ,'(3x,A40,I16)')adjl('Number of Dirichlet (q=0) faces:',40), n_dirichlet_faces
    write(6  ,'(3x,A39)',advance='no')'---face ids:                               '
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
    write(iow,'(3x,A40)')adjl('There are no Dirichlet (q=0) faces.',40)
    write(6  ,'(3x,A40)')adjl('There are no Dirichlet (q=0) faces.',40)
endif

if (log_convergence_scheme) then
    if (scheme_type.eq.1 .or. scheme_type.eq.2) then
        write(iow,'(3x,A40,I16)')adjl('Convergence scheme chosen:',40), scheme_type
        write(6  ,'(3x,A40,I16)')adjl('Convergence scheme chosen:',40), scheme_type
    else
        write(ERROR_MESSAGE,'(''Convergence scheme does not exist! Please choose a value between 1 and 2'',I10)') scheme_type
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    scheme_type = 2
    write(iow,'(3x,A40)')adjl('convergence scheme not found..',40), scheme_type
    write(iow,'(3x,A40,I16)')adjl('---It was set to the default value:',40), scheme_type
    write(6  ,'(3x,A40)')adjl('convergence scheme not found..',40), scheme_type
    write(6  ,'(3x,A40,I16)')adjl('---It was set to the default value:',40), scheme_type
endif


if (log_mumps_matrix_type) then
    if (mumps_matrix_type == 0) then
        write(iow,'(3x,A40,1x,A14,I1,A1)')adjl('mumps matrix type:',40),'nonsymmetric (', mumps_matrix_type,')'
        write(6  ,'(3x,A40,1x,A14,I1,A1)')adjl('mumps matrix type:',40),'nonsymmetric (', mumps_matrix_type,')'
    elseif (mumps_matrix_type == 1) then
        write(iow,'(3x,A40,1x,A21,I1,A1)')adjl('mumps matrix type:',40),'symmetric pos. def. (', mumps_matrix_type,')'
        write(6  ,'(3x,A40,1x,A21,I1,A1)')adjl('mumps matrix type:',40),'symmetric pos. def. (', mumps_matrix_type,')'
    elseif (mumps_matrix_type == 2) then
        write(iow,'(3x,A40,1x,A18,I1,A1)')adjl('mumps matrix type:',40),'general symmetric(', mumps_matrix_type,')'
        write(6  ,'(3x,A40,1x,A18,I1,A1)')adjl('mumps matrix type:',40),'general symmetric(', mumps_matrix_type,')'
    else
        write(ERROR_MESSAGE,'(''Incorrect mumps matrix type..'',I16)') mumps_matrix_type
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
else
    mumps_matrix_type = 0
    write(iow,'(3x,A40)')adjl('mumps matrix type not detected..',40)
    write(iow,'(3x,A40,I10)')adjl('---It was set to the nonsymmetric:',40),mumps_matrix_type
    write(6  ,'(3x,A40)')adjl('mumps matrix type not detected..',40)
    write(6  ,'(3x,A40,I10)')adjl('---It was set to the nonsymmetric:',40),mumps_matrix_type
endif

if (log_read_field.and.readfield==1) then
    write(iow,'(A43,A15)')adjl('*Field will be read from file:',40),field_filename
    write(6  ,'(A43,A15)')adjl('*Field will be read from file:',40),field_filename
endif
if (.not.log_read_field.or.readfield.ne.1) then
    write(iow,'(/A40)')adjl('*Field will be initialized to zero..',40)
    write(6  ,'(/A40)')adjl('*Field will be initialized to zero..',40)
endif
!*******************************************************************!
!                Initialize some useful quantinties                 !
!*******************************************************************!

! Hamaker constants
Aps   = Aps*1.e-20
Asio2 = Asio2*1.e-20

write(iow,'(/''*Initialization of usefull quantities..'')')
write(6  ,'(/''*Initialization of usefull quantities..'')')
ds = 1.d0/dble(ns)
write(iow,'(3x,A40,E16.9)')adjl('ds:',40),ds
write(6  ,'(3x,A40,E16.9)')adjl('ds:',40),ds

! Calculate the radius of gyration
Rgyr = 1.54d00 * dsqrt(CN * (chainlen)/6.d00)
write(iow,'(3x,A40,E16.9,A10)')adjl('Radius of Gyration:',40),Rgyr,' Angstroms'
write(6  ,'(3x,A40,E16.9,A10)')adjl('Radius of Gyration:',40),Rgyr,' Angstroms'


!calculate segment density in the bulk rho_0 in mol_segments/m3
rho_0 = chainlen*massden/(chainlen*mon_mass )*1.d06

write(iow,'(3x,A40,E16.9,A8)')adjl('Segment density in bulk (rho):',40),rho_0,' mol/m^3'
write(6  ,'(3x,A40,E16.9,A8)')adjl('Segment density in bulk (rho):',40),rho_0,' mol/m^3'

kapa = chainlen/(kappa_T * boltz_const_Joule_molK * Temp * rho_0)

write(iow,'(3x,A40,E16.9)')adjl('kapa = 1/[k_T k_B T rho_0]:',40),kapa
write(iow,'(3x,A40,E16.9)')adjl('kapa / chainlen:',40),kapa / chainlen
write(6  ,'(3x,A40,E16.9)')adjl('kapa = 1/[k_T k_B T rho_0]:',40),kapa
write(6  ,'(3x,A40,E16.9)')adjl('kapa / chainlen:',40),kapa / chainlen

return

!--------------------------------------------------------------------------------!
end subroutine scfinout
