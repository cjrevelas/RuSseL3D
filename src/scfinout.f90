subroutine scfinout
!--------------------------------------------------------------------------------!   
use xdata
use mdata
!--------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------!    
open(unit=iow, file = 'scft_out.txt')   
open(unit=ior, file = 'gaussdat.txt')
 
write(iow,'(''Polymer FEM Bulk v.1.2  27 Jun 19 '')')
write(iow,'(''numnp='',I6,'' ns='',I6 )') numnp, ns
write(iow,'(''Polumer Bulk with solid surfaces   ''/  &
            ''---------------------------------''/    &
            ''SCF theory, Gaussian string model''/    &
            ''Finite Element solution with     ''/    &
            ''successive substitutions         '')')
!***********************************************************!
!                      CHOOSE THE MODEL                     !
!***********************************************************!
!write(6,*) '0-> Helfand Model'
!write(6,*) 'Choose gaussdat or Sanchez Lacombe:'
!write(6,*) '0-> gaussdat params'
!read (6,*) slparams
slparams=0

!Read seed for random number generator NOT USED
read(ior,'(I9)') iseed
write(iow,'(/A30,I11)') '1->Sphere 0-> Film', iseed

write(iow,'(/A30,I11)') 'Number of nodes ',    numnp
write(iow,'(/A30,I11)') 'Number of elements ', numel

read(ior,'(I6)') ns
write(iow,'(/A30,I11)') 'Contour length (s) discretization ', ns

!Check to make sure ns is even
if (mod(ns,2).ne.0) then
    write(6,'(''ns not even: '',I6)') ns
    stop
endif

allocate(wa(numnp),wa_new(numnp),wa_in(numnp),Ufield(numnp))

allocate(qf(numnp,2),qf_final(numnp,ns+1))
      
allocate(phia_new(numnp),phi(numnp))
      
!Interval used for dis!retization of dimensionless contour length variable s
ds = 1.d0/dble(ns)

!Read in chain length (in monomer units)
read(ior,'(I10)') chainlen
write(iow,'(/A30,I11)') 'Chain length ', chainlen

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
Rgyr = 1.54d00 * dsqrt(CN * (dfloat(chainlen))/6.d00)

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

!Read additional parameters
read(ior,'(E16.9)') mon_mass
write(iow,'(/A30,F16.4)')'Mass  of a monomer', mon_mass
 
read(ior,'(E16.9)') sphere_radius
write(iow,'(/A30,F16.4)' )' sphere_radius', sphere_radius
 
read(ior,'(E16.9)') sigma1
write(iow,'(/A30,F16.4,    '' Angstroms'')') adjustl('sigma polymer'), sigma1
 
read(ior,'(E16.9)')sigma2
write(iow,'(/A30,F16.4,    '' Angstroms'')')adjustl('sigma solid'), sigma2
 
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
      
!calculate segment density in the bulk rho_0 in mol_segments/m3
rho_0 = dfloat(chainlen)*massden/(dfloat(chainlen)*mon_mass )*1.d06
    
write(6,*) 'rho_0 (mol/m^3)'
write(6,*) rho_0
write(iow,'(/A30,F16.4,'' mol/m3'')') 'Segment density in bulk', rho_0

kapa = dfloat(chainlen)/(kappa_T * k_B * Temp * rho_0)
write(6,*) 'kapa'
write(6,*) kapa

write(iow,'(/A30,F16.4)') 'kapa = 1/[k_T k_B T rho_0]', kapa/ dfloat(chainlen)

close(ior)

return
!--------------------------------------------------------------------------------!
end subroutine scfinout
       
       
       
       
