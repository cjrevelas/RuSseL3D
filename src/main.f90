program FEM_3D
    
    
    use xdata
    use mdata
    use kcw
  

    implicit none

    real*8 surf_pot

!===================================================================!
!
!             Gaussdat Parameter values
!
!===================================================================!
         ior=55
         iow=10    
         
         call mesh_io_3d
         
         
         print*,'numnp',numnp,'numel',numel
       
         call scfinout  

!-------------------------------------------------------------------!
!      Initialize Simpson coeffi!ients
!     for Contour integration
!-------------------------------------------------------------------!

         call simpsonkoef_s
!checked

         open(unit=211, file = 'Usolid_out.txt')

	     do k1 = 1, numnp
		      distance =xc(1,k1)
              Ufield(k1) =surf_pot(distance)
              write(211,('(e16.9,2x,e19.9)'))distance ,Ufield(k1) 
         end do!k1

!-------------------------------------------------------------------!
!
!                       Initialize fields
!
!-------------------------------------------------------------------!
!
!........This overwrites the fields with zeros

         do i1 =1, numnp
              wa(i1) = 0.d0
         end do
!........Field values read in
         if (show.eq.1) then
              open(unit = 21, file = 'field_in.txt')
              do i1 =1, numnp
                   read (21,'(18X,E16.9)') wa(i1)
                   wa(i1)=wa(i1)!*real(chainlen)
              end do
              close(21)
         end if

!-------------------------------------------------------------------!
!
!                       LOOPS FOR SOLUTION                          !
!
!-------------------------------------------------------------------!

         write(*,*) 'Iteration,  Adh. tension (mN/m),  error(beta N w)'
         write(iow,'(a10,2X,a16,2X,a16,2X,a16)')'kk','adh_ten_alt','error'

         kk=0
         error=200000.
         do while ((kk.lt.iterations).and.(error.gt.max_error))
              kk=kk+1
        
         call Matrix_assemble(c,k,w)  
        print*,'edw kaleite  h edwards'
         call edwards_free_film_fem 


              if (mod(kk,10000).eq.0) then  
                   open(unit=21, file = 'field_out.txt')
                   do k1 =1,numnp
                        write(21,'(E16.9,2X,E16.9)') xc(1,k1), wa(k1)
                   end do
                   close(21)
              end if
                   
         call part_fun_phi                  
        
         
        
!-------------------------------------------------------------------!
!..................Periodic output of Profiles
!-------------------------------------------------------------------!
              open (file='reduced_density_profiles.txt', unit=120)
              if (mod(kk,1000).eq.0.d0) then
                   write(120,'(a16,2X,a16)') 'z','phi(z)'
                   do k1 = 1,numnp
                         write (120,'(E16.9,2X,E16.9)') xc(1,k1), phia_new(k1)
                   end do
              end if
             close(120) 
           
         call adhesion_tension_alternative

!-------------------------------------------------------------------!
!*******************************************************************!
!***********************New field calculations**********************!
!*******************************************************************!
!*******************************************************************!

		      do k1 =1,numnp
			       distance = xc(1,k1)
	               wa_new(k1) = kapa * (phia_new(k1)- 1.d0) + Ufield(k1)
               end do
!-------------------------------------------------------------------!
!........Test for convergence on fields, using absolute error
!-------------------------------------------------------------------!
!........Norm of deviation between fields in two su!!essive
!........iterations accumulated in error                            !
!
		      error = 0.d00
              do k1 = 1,numnp
                   error = max(error, dabs(wa_new(k1) - wa(k1)))
              end do!k1
!
!............Mix in fraction of new field
!............For field expression, see eq (8a) of
!............Daoulas et al. 2005.		
              do k1 =1,numnp 
                   wa(k1) = (1.d0-fraction) * wa(k1) + fraction* wa_new(k1)
              end do!k1
!
              if (mod(kk,1000).eq.0d0)then
                   open(123,file='field_compare.txt')
                   do k1=1,numnp
                         write(123,'(i10,2x,e16.9,2x,e16.9)')k1,wa(k1),wa_new(k1)
                   end do
		       end if
!
              if (mod(kk,100).eq.0.d0) then
		            write(*,*)   kk, adh_ten_alt, error
	                write(iow,'(I10,2X,E16.9,2X,E16.9)') kk, adh_ten_alt,error
              end if

         end do!kk
     

         if (pr_on==1)then 
              call qprint
         end if


	     if (error.lt.max_error)then
              write(iow,'(/''Convergence of max error'',f16.9)')error
         else
              write(iow,'(/''Convergence of '',i10, '' iterations'')')iterations
         end if
         
	     write(iow,'(''-----------------------------------'')')
	     write(iow,'(''Adhesion tension (mN/m) '',E16.9)') adh_ten_alt
	     write(iow,'(''Partition function Q =  '',E16.9)')part_func
         write(iow,'(''            n/n_bulk =  '',E16.9)') nch_per_area *dfloat(chainlen) / (rho_0 * lx  * 1.d-10)

end
      
