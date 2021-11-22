subroutine scfinout
    
    use xdata
    use mdata

    implicit none
    
         
         open(unit=iow, file = 'scft_out.txt')   
         open(unit=ior, file = 'gaussdat.txt')
      
!
         write(iow,'(''Polymer FEM Bulk v.1.2  27 Jun 19 '')')
         write(iow,'(''maxnx='',i6,'' maxxns='',i6 )') maxnx,maxns
         write(iow,'(''Polumer Bulk with solid surface   ''/         &
                      ''---------------------------------''/         &
                      ''SCF theory, Gaussian string model''/         &
                      ''Finite Element solution with     ''/         &
                      ''successive substitutions         '')')
!=======================================================================
!
!            !hose the Model
!
!=======================================================================
!!hoose your model San!ez La!ombe model or Helfand

         write(*,*) 'Grafted polymer Bulk v.1.2  01 Jun 19'
!    write(*,*) 'Choose San!ez-Lacombe or Helfand :'

         write(*,*) '0-> Helfand Model'
!

!    write(*,*) '!hoose gaussdat or San!hez La!ombe  :'

         write(*,*) '0-> gaussdat params'
!    read (*,*) slparams
         slparams=0


    
    
    
    !Read seed for random number generator NOT USED
	     read(ior,'(I9)') iseed
         write(iow,'(/a30,I11)') '1->Sphere 0-> Film',iseed
        
!
!Domain of height lx in the dire!tion normal to the
!interfa!e is !onsidered.
!Read in lx, in units of Angstroms.
!Also read in nx, the number of intervals to be used
!in the dis!retization of lx.
!Note: nx must be an even integer.
      
         write(iow,'(/a30,f16.4)') 'Domain  in Angstroms ',lx
         write(iow,'(/a30,I11)') 'Number of nodes ', numnp
         write(iow,'(/a30,I11)') 'Number of elements ', numel

!

!
!Read in ns, the number of intervals to be used
!in the discretization of the dimensionless contour.
!Note: ns must be an even number.
	     read(ior,'(I6)') ns
         write(iow,'(/a30,I11)') 'Contour length (s) discretization ', ns
      
!
!!he!k to make sure ns is even
	     if(mod(ns,2).ne.0) then
		      write(*,'(''ns not even: '',I6)') ns
              stop
         end if
!
!Interval used for dis!retization of dimensionless
!!ontour length variable s
	     ds = 1.d0/dble(ns)
!
!Read in !hain length (in monomer units, i.e. methylenes/methyls)
	     read(ior,'(I10)') chainlen
         write(iow,'(/a30,I11)') 'Chain length ',chainlen
!
!
!
!Read absolute temperature, in K.
	     read(ior,'(E16.9)') Temp
         write(iow,'(/a30,f16.4,'' K'')') 'Temperature ',Temp
!
!Read redu!ed well depth of surfa!e potential, Vored,
!in units of k_B T.
!Surfa!e potential assumed to be of a square well form.
         read(ior,'(E16.9)') Aps
         write(iow,'(/a30,f16.4,    '' e-20 Joule'')')adjustl('Hamaker constant of polymer'), Aps
      
	     read(ior,'(E16.9)') Asio2
         write(iow,'(/a30,f16.4,    '' e-20 Joule'')')adjustl('Hamaker constant of solid') , Asio2
       
         Aps= Aps*1.e-20  
      
         Asio2=Asio2*1.e-20   
    
                

!Read mass density in the bulk, in g/!m3
	     read(ior,'(E16.9)') massden
         write(iow,'(/a30,f16.4,'' g/cm3'')') 'Mass density', massden

!Read isothermal !ompressibility kappa_T, in Pa^-1
         read(ior,'(E16.9)') kappa_T
         write(iow,'(/a30,f16.4,   ''x10-9 Pa^-1'')') 'Isothermal compressibility',kappa_T*1.e9
            
!Read Characteristic ratio CN
         read(ior,'(E16.9)') CN
         write(iow,'(/a30,f16.4)') 'Chain characteristic ratio C_N ' ,CN

	     Rgyr = 1.54d00 * dsqrt(CN * (dfloat(chainlen))/6.d00)
!
         write(*,*) 'The gyration radius, in Angstroms, is:'
         write(*,*)  Rgyr
         write(iow,'(/a30,f16.4,    '' Angstroms'')') 'radius of gyration R_g' ,Rgyr
         
!
!
!
!Read maximum number of iterations and toleran!e for !onvergen!e
         read(ior,'(I10,E16.9)') iterations, max_error
         write(iow,'(/a30,I11)') 'Maximum iterations ',iterations
         write(iow,'(/a30,f16.4)')'Tolerance for convergence of field ',max_error

!Read fra!tion of new field to be mixed in with old field
!in damped iterations
	     read(ior,'(E16.9)') fraction
         write(iow,'(/a30,f16.4)'  ) 'Fraction of new field mixed in with old field in damped iterations ',fraction
         read(ior,'(E16.9)') mon_mass
         write(iow,'(/a30,f16.4)')'Mass  of a monomer', mon_mass
         read(ior,'(E16.9)') sphere_radius
         write(iow,'(/a30,f16.4)' )' sphere_radius',sphere_radius
         read(ior,'(E16.9)') sigma1
         write(iow,'(/a30,f16.4,    '' Angstroms'')') adjustl('sigma polymer'),sigma1
         read(ior,'(E16.9)')sigma2
         write(iow,'(/a30,f16.4,    '' Angstroms'')')adjustl('sigma solid'), sigma2
         read(ior,'(i10)',iostat=lshow ) show
        

         if (lshow <0)then 
               show=0
               write(iow,'(/a30,a16)')adjustl('Initial field'), 'DEFAULT FIELDIN'
         else 
               if (show==1)then
                  write(iow,'(/a30,a16)')adjustl('Initial field'), 'FIELD_IN'
               else
                   write(iow,'(/a30,a16)')adjustl('Initial field'), 'FIELD ZERO'
               end if 
         end if 
        read(ior,'(i10)',iostat=lshow ) pr_on
        
!       print*, pr_on
         if (lshow <0)then 
              pr_on=0
              write(iow,'(/a30,a16)')adjustl('QPRINT'), 'DEFAULT OFF'
         else 
              if (pr_on==1)then
                   write(iow,'(/a30,a16)')adjustl('QPRINT'), 'ON'
              else 
                   write(iow,'(/a30,a16)')adjustl('QPRINT'), 'OFF'
              end if 
         end if 
      
!!al!ulate segment density in the bulk, rho_0,
!in mol segments per !ubi!meter
         rho_0 = dfloat(chainlen)*massden/(dfloat(chainlen)*mon_mass )*1.d06
!
     
         write(*,*) 'rho_0 (mol/m^3)'
         write(*,*) rho_0
         write(iow,'(/a30,f16.4,'' mol/m3'')')'Segment density in bulk', rho_0
        
!The !ode is made to be !ompatible with previous gauss11
!The !ompressibility is only for Helfand Model
!kapa is referred  only in Helfand Model
	     kapa = dfloat(chainlen)/(kappa_T * k_B * Temp * rho_0)
!
         write(*,*) 'kapa'
         write(*,*) kapa


         write(iow,'(/a30,f16.4)') 'kapa = 1/[k_T k_B T rho_0] ',kapa/ dfloat(chainlen)
     
! FOR FD FORMULATION     
!      diff_number(1)  = ds*Rgyr*Rgyr/(2.d0*dx*dx)
!      write(*,*) 'diff_number'
!      write(*,*) diff_number(1)
!      write(iow,'(/a30,f16.4)') 'Dimensionless diffusion number ',diff_number(1)

         close (ior)

    return

end 
       
       
       
       