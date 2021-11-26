subroutine scfinout
!--------------------------------------------------------------------------------!   
use xdata
use constants
use mdata
!--------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------!    
open(unit=iow, file = 'scft.out.txt')   
open(unit=ior, file = 'gaussdat.in.txt')

write(iow,'(''Polymer FEM Bulk v.1.2  27 Jun 19 '')')
write(iow,'(''numnp='',I6,'' ns='',I6 )') numnp, ns
write(iow,'(''Polumer Bulk with solid surfaces   ''/  &
          & ''---------------------------------''/    &
          & ''SCF theory, Gaussian string model''/    &
          & ''Finite Element solution with     ''/    &
          & ''successive substitutions         '')')
!***********************************************************!
!                      CHOOSE THE MODEL                     !
!***********************************************************!

!Read seed for random number generator NOT USED
read(ior,'(I9)') iseed
write(iow,'(/A30,I11)') '1->Sphere 0-> Film', iseed

write(iow,'(/A30,I11)') 'Number of nodes ',    numnp
write(iow,'(/A30,I11)') 'Number of elements ', numel

read(ior,'(I6)') ns
write(iow,'(/A30,I11)') 'Contour length (s) discretization ', ns

!Check to make sure ns is even
if (mod(ns,2).ne.0) then
    write(6,'(''ns is not an even number: '',I6)') ns
    stop
endif

allocate(wa(numnp),wa_new(numnp),wa_mix(numnp),Ufield(numnp))

allocate(qf(numnp,2),qf_final(numnp,ns+1))

allocate(phia_new(numnp),phi(numnp))

!Interval used for discretization of dimensionless contour length variable s
ds = 1.d0/dble(ns)

!Read in chain length (in monomer units)
read(ior,'(E16.9)') chainlen
write(iow,'(/A30,E16.9)') 'Chain length ', chainlen

!Read absolute temperature in K.
read(ior,'(E16.9)') Temp
write(iow,'(/A30,F16.4,'' K'')') 'Temperature ', Temp

!Read Hamaker potential parameters
read(ior,'(E16.9)') Aps
write(iow,'(/A30,F16.4,    '' e-20 Joule'')')adjustl('Hamaker constant of polymer'), Aps

read(ior,'(E16.9)') Asio2
write(iow,'(/A30,F16.4,    '' e-20 Joule'')')adjustl('Hamaker constant of solid') , Asio2

Aps   = Aps*1.e-20
Asio2 = Asio2*1.e-20

!Read mass density in the bulk in g/cm3
read(ior,'(E16.9)') massden
write(iow,'(/A30,F16.4,'' g/cm3'')') 'Mass density', massden

!Read isothermal compressibility kappa_T in Pa^-1
read(ior,'(E16.9)') kappa_T
write(iow,'(/A30,F16.4,   ''x10-9 Pa^-1'')') 'Isothermal compressibility', kappa_T*1.e9
   
!Read Characteristic ratio CN
read(ior,'(E16.9)') CN
write(iow,'(/A30,F16.4)') 'Chain characteristic ratio C_N ', CN

!Calculate the radius of gyration
Rgyr = 1.54d00 * dsqrt(CN * (chainlen)/6.d00)

write(6,*) 'The gyration radius, in Angstroms, is:'
write(6,*)  Rgyr
write(iow,'(/A30,F16.4,    '' Angstroms'')') 'radius of gyration R_g', Rgyr

!Read maximum number of iterations and tolerance for convergence
read(ior,'(I10,E16.9)') iterations, max_error
write(iow,'(/a30,I11)') 'Maximum iterations ', iterations
write(iow,'(/A30,F16.4)')'Tolerance for convergence of field ', max_error

!Read fraction of new field to be mixed in with old field in damped iterations
read(ior,'(E16.9)') fraction
write(iow,'(/a30,f16.4)'  ) 'Fraction of new field mixed in with old field in damped iterations ',fraction

if (fraction > 1.d0) then
    write(iow,'(/a30,f16.4)'  ) 'ERROR: fraction is larger than one!'
    STOP
endif

!Read additional parameters
read(ior,'(E16.9)') mon_mass
write(iow,'(/A30,F16.4)')'Mass  of a monomer', mon_mass

read(ior,'(E16.9)') sphere_radius
write(iow,'(/A30,F16.4)' )' sphere_radius', sphere_radius

read(ior,'(E16.9)') sigma1
write(iow,'(/A30,F16.4,    '' Angstroms'')') adjustl('sigma polymer'), sigma1

read(ior,'(E16.9)') sigma2
write(iow,'(/A30,F16.4,    '' Angstroms'')')adjustl('sigma solid'), sigma2

!Read field
read(ior,'(I10)',iostat=lshow) show
if (lshow <0)then
    show=0
    write(iow,'(/A30,A16)')adjustl('Initial field'), 'DEFAULT FIELDIN'
else
    if (show==1)then
       write(iow,'(/A30,A16)')adjustl('Initial field'), 'FIELD_IN'
    else
       write(iow,'(/A30,A16)')adjustl('Initial field'), 'FIELD ZERO'
    endif
endif

!print q
read(ior,'(I10)',iostat=lshow ) pr_on

if (lshow<0)then
    pr_on = 0
    write(iow,'(/A30,A16)')adjustl('QPRINT'), 'DEFAULT OFF'
else
    if (pr_on==1) then
        write(iow,'(/A30,A16)')adjustl('QPRINT'), 'ON'
    else
        write(iow,'(/A30,A16)')adjustl('QPRINT'), 'OFF'
    endif
endif

!read mix coeff frac
read(ior,'(E16.9)',iostat=lshow ) mix_coef_frac

if (lshow<0)then
    mix_coef_frac = 0
    write(iow,'(/A30,E16.9)')adjustl('mix_coef_frac'), mix_coef_frac
else
    write(iow,'(/A30,E16.9)')adjustl('mix_coef_frac'), mix_coef_frac
endif


!read mix coeff kapa
read(ior,'(E16.9)',iostat=lshow ) mix_coef_kapa

if (lshow<0)then
    mix_coef_frac = 1
    write(iow,'(/A30,E16.9)')adjustl('mix_coef_kapa'), mix_coef_kapa
else
    write(iow,'(/A30,E16.9)')adjustl('mix_coef_kapa'), mix_coef_kapa
endif

write(*,*)"KAPAAA",mix_coef_kapa

!read dirichlet faces
read(ior,'(I10)',iostat=lshow ) n_dirichlet_faces

if (lshow<0)then
    n_dirichlet_faces = 0
    write(iow,'(/A30,I10)')adjustl('number faces applying Dirichlet BC, q=0?'), n_dirichlet_faces
else
    write(iow,'(/A30,I10)')adjustl('number faces applying Dirichlet BC, q=0?'), n_dirichlet_faces
    allocate(ids_dirichlet_faces(n_dirichlet_faces))
endif

if (n_dirichlet_faces>0) then
    write(iow,'(/A30,I10)',advance='no')adjustl('ids..')
    do i1 = 1, n_dirichlet_faces
        read(ior,'(I3)',advance='no') id
        write(iow,'(I3)',advance='no') id
        ids_dirichlet_faces(i1) = id
    enddo
    read(ior,*)
    write(iow,*)
endif

!choose convergence scheme
read(ior,'(I10)',iostat=lshow ) scheme_type

if (lshow<0)then
    scheme_type = 0
    write(iow,'(/A30,I10)')adjustl('scheme_type'), scheme_type
else
    write(iow,'(/A30,I10)')adjustl('scheme_type'), scheme_type
endif

!STOP

!calculate segment density in the bulk rho_0 in mol_segments/m3
rho_0 = chainlen*massden/(chainlen*mon_mass )*1.d06

write(6,*) 'rho_0 (mol/m^3)'
write(6,*) rho_0
write(iow,'(/A30,F16.4,'' mol/m3'')') 'Segment density in bulk', rho_0

kapa = chainlen/(kappa_T * boltz_const_Joule_molK * Temp * rho_0)
write(6,*) 'kapa'
write(6,*) kapa

write(iow,'(/A30,F16.4)') 'kapa = 1/[k_T k_B T rho_0]', kapa / chainlen

close(ior)

return
!--------------------------------------------------------------------------------!
end subroutine scfinout
