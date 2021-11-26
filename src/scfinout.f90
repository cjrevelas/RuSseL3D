subroutine scfinout
!--------------------------------------------------------------------------------!   
use xdata
use constants
use mdata
!--------------------------------------------------------------------------------!
implicit none

CHARACTER(100) :: line
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

open(unit=iow, file = 'scft.out.txt')   
open(unit=ior, file = 'gaussdat.in.txt')


write(iow,'(''Polymer FEM Bulk v.1.2  27 Jun 19 '')')
write(iow,'(''numnp='',I6,'' ns='',I6 )') numnp, ns
write(iow,'(''Polumer Bulk with solid surfaces   ''/  &
          & ''---------------------------------''/    &
          & ''SCF theory, Gaussian string model''/    &
          & ''Finite Element solution with     ''/    &
          & ''successive substitutions         '')')


do
    read(ior,'(A100)',IOSTAT=Reason) line

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
            read(line,'(F16.4)') max_error
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
            read(line,'(I10)') show
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
        elseif (index(line,'# n dirichlet faces') > 0) then
            read(line,'(I10)') n_dirichlet_faces
            allocate(ids_dirichlet_faces(n_dirichlet_faces))
            do i1 = 1, n_dirichlet_faces
                 read(ior,'(I10)')id
                 ids_dirichlet_faces(i1) = id
            enddo
            log_n_dirichlet_faces = .true.
        endif
   endif
enddo

close(ior)

write(iow,'(/A30,I11)') '1->Sphere 0-> Film', iseed
write(iow,'(/A30,I11)') 'Number of nodes ',    numnp
write(iow,'(/A30,I11)') 'Number of elements ', numel
write(iow,'(/A30,I11)') 'Contour length (s) discretization ', ns
if (mod(ns,2).ne.0) then
    write(6,'(''ns is not an even number: '',I6)') ns
    stop
endif
ds = 1.d0/dble(ns)
!--------------------------------------------------------------------------------!    
allocate(wa(numnp),wa_new(numnp),wa_mix(numnp),Ufield(numnp))
allocate(qf(numnp,2),qf_final(numnp,ns+1))
allocate(phia_new(numnp),phi(numnp))
write(iow,'(/A30,E16.9)') 'Chain length ', chainlen
write(iow,'(/A30,F16.4,'' K'')') 'Temperature ', Temp
write(iow,'(/A30,F16.4,    '' e-20 Joule'')')adjustl('Hamaker constant of polymer'), Aps
write(iow,'(/A30,F16.4,    '' e-20 Joule'')')adjustl('Hamaker constant of solid') , Asio2
Aps   = Aps*1.e-20
Asio2 = Asio2*1.e-20
write(iow,'(/A30,F16.4,'' g/cm3'')') 'Mass density', massden

write(iow,'(/A30,F16.4,   ''x10-9 Pa^-1'')') 'Isothermal compressibility', kappa_T*1.e9
write(iow,'(/A30,F16.4)') 'Chain characteristic ratio C_N ', CN
!
!!Calculate the radius of gyration
Rgyr = 1.54d00 * dsqrt(CN * (chainlen)/6.d00)
!
write(6,*) 'The gyration radius, in Angstroms, is:'
write(6,*)  Rgyr
write(iow,'(/A30,F16.4,    '' Angstroms'')') 'radius of gyration R_g', Rgyr
!
!!Read maximum number of iterations and tolerance for convergence
!read(ior,'(I10,E16.9)') iterations, max_error
write(iow,'(/a30,I11)') 'Maximum iterations ', iterations
write(iow,'(/A30,F16.4)')'Tolerance for convergence of field ', max_error
!
!!Read fraction of new field to be mixed in with old field in damped iterations
!read(ior,'(E16.9)') fraction
write(iow,'(/a30,f16.4)'  ) 'Fraction of new field mixed in with old field in damped iterations ',fraction
!
if (fraction > 1.d0) then
    write(iow,'(/a30,f16.4)'  ) 'ERROR: fraction is larger than one!'
    STOP
endif
!
!!Read additional parameters
write(iow,'(/A30,F16.4)')'Mass  of a monomer', mon_mass
!
!read(ior,'(E16.9)') sphere_radius
write(iow,'(/A30,F16.4)' )' sphere_radius', sphere_radius
!
!read(ior,'(E16.9)') sigma1
write(iow,'(/A30,F16.4,    '' Angstroms'')') adjustl('sigma polymer'), sigma1
!
!read(ior,'(E16.9)') sigma2
write(iow,'(/A30,F16.4,    '' Angstroms'')')adjustl('sigma solid'), sigma2
!
!Read field
if (log_read_field) then
    if (show==1)then
       write(iow,'(/A30,A16)')adjustl('Initial field'), 'FIELD_IN'
    else
       write(iow,'(/A30,A16)')adjustl('Initial field'), 'FIELD ZERO'
    endif
else
    show=0
    write(iow,'(/A30,A16)')adjustl('Initial field'), 'DEFAULT FIELDIN'
endif


write(iow,'(/A30,A16)')adjustl('QPRINT'), 'ON'

if (log_mix_coef_fraction) then
    write(iow,'(/A30,E16.9)')adjustl('mix_coef_frac'), mix_coef_frac
else
    mix_coef_frac = 0
    write(iow,'(/A30,E16.9)')adjustl('mix_coef_frac'), mix_coef_frac
endif

if (log_mix_coef_kapa) then
    write(iow,'(/A30,E16.9)')adjustl('mix_coef_kapa'), mix_coef_kapa
else
    mix_coef_frac = 1
    write(iow,'(/A30,E16.9)')adjustl('mix_coef_kapa'), mix_coef_kapa
endif

if (log_n_dirichlet_faces) then
    write(iow,'(/A30,I10)')adjustl('number faces applying Dirichlet BC, q=0?'), n_dirichlet_faces
    write(iow,'(/A30,I10)',advance='no')adjustl('ids..')
    do i1 = 1, n_dirichlet_faces
        write(iow,'(I3)',advance='no') ids_dirichlet_faces(i1)
    enddo
    write(iow,*)
else
    n_dirichlet_faces = 0
    allocate(ids_dirichlet_faces(1))
    ids_dirichlet_faces(1)=-1
    write(iow,'(/A30,I10)')adjustl('number faces applying Dirichlet BC, q=0?'), n_dirichlet_faces
endif

!
!!choose convergence scheme

if (log_convergence_scheme) then
    write(iow,'(/A30,I10)')adjustl('scheme_type'), scheme_type
else
    scheme_type = 1
    write(iow,'(/A30,I10)')adjustl('scheme_type'), scheme_type
endif

!calculate segment density in the bulk rho_0 in mol_segments/m3
rho_0 = chainlen*massden/(chainlen*mon_mass )*1.d06

write(6,*) 'rho_0 (mol/m^3)'
write(6,*) rho_0
write(iow,'(/A30,F16.4,'' mol/m3'')') 'Segment density in bulk', rho_0

kapa = chainlen/(kappa_T * boltz_const_Joule_molK * Temp * rho_0)
write(6,*) 'kapa'
write(6,*) kapa

write(iow,'(/A30,F16.4)') 'kapa = 1/[k_T k_B T rho_0]', kapa / chainlen

!STOP

return
!--------------------------------------------------------------------------------!
end subroutine scfinout
