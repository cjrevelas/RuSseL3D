program FEM_3D
!------------------------------------------------------------------------------------------------------------------!     
use xdata
use mdata
use kcw
!------------------------------------------------------------------------------------------------------------------!  
implicit none
!------------------------------------------------------------------------------------------------------------------!
real(8) :: surf_pot
!------------------------------------------------------------------------------------------------------------------!
!*******************************************************************!
!                   GAUSSDAT PARAMETER VALUES                       !
!*******************************************************************!
ior = 55
iow = 10
   
open(unit=54, file='times.txt')

call CPU_TIME(start)  
call mesh_io_3d
call CPU_TIME(t1)

write(54, '("Time meshio = ",I6.3,"minutes",F6.3," seconds.")') int(t1-start)/60, mod((t1-start),60.)

write(54,*) 'numnp', numnp, 'numel', numel
write(6,*)  'numnp', numnp, 'numel', numel
 
call scfinout  
!*******************************************************************!
!                 INITIALIZE SIMPSON COEFFICIENTS                   !
!                      FOR CONTOUR INTEGRATION                      !
!*******************************************************************!
call simpsonkoef_s

open(unit=211, file = 'Usolid_out.txt')

do k1 = 1, numnp
     distance   = xc(1,k1)
     Ufield(k1) = surf_pot(distance)
     write(211,('(E16.9,2X,E19.9)')) distance, Ufield(k1) 
enddo
!*******************************************************************!
!                       INITIALIZE FIELDS                           !
!*******************************************************************!
do i1 = 1, numnp
     wa(i1) = 0.d0
enddo

if (show.eq.1) then
   open(unit=21, file = 'field_in.txt')
   do i1 = 1, numnp
        read (21,'(18X,E16.9)') wa(i1)
        wa(i1)=wa(i1) !*real(chainlen)
   enddo
   close(21)
endif
!*******************************************************************!
!                       LOOPS FOR SOLUTION                          !
!*******************************************************************!
write(6,*) 'Iteration,  Adh. tension (mN/m),  error(beta N w)'
write(iow,'(A10,2X,A16,2X,A16,2X,A16)') 'kk', 'adh_ten_alt', 'error'

kk=0
error=200000.

do while ((kk.lt.iterations).and.(error.gt.max_error))
    kk=kk+1

    call CPU_TIME(t2)
    call matrix_assemble   
    call CPU_TIME(t3)

    write(54, '("Time of matrix assembly = ",I6.3,"minutes",F6.3," seconds.")')  int(t3-t2)/60, mod((t3-t2),60.)
  
    call edwards_free_film_fem 
    call CPU_TIME(t4)

    write(54, '("Time of Edwards solution = ",I6.3,"minutes",F6.3," seconds.")') int(t4-t3)/60, mod((t4-t3),60.)
   
    write(54, '(e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,e16.9)') start, t1, t2, t3, t4
   
    close(54)
   
    if (mod(kk,10000).eq.0) then  
        open(unit=21, file = 'field_out.txt')
        do k1 = 1, numnp
            write(21,'(E16.9,2X,E16.9)') xc(1,k1), wa(k1)
        enddo
        close(21)
    endif
  
    call part_fun_phi                  
    !*******************************************************************!
    !                     PERIODIC EXPORT OF PROFILES                   !
    !*******************************************************************!
    open (unit=120, file = 'reduced_density_profiles.txt')
    !if (mod(kk,1000).eq.0.d0) then
        write(120,'(A16,2X,A16)') 'z','phi(z)'
        do k1 = 1, numnp
            write(120,'(E16.9,2X,E16.9,2X,E16.9,2X,E16.9)') xc(1,k1), xc(2,k1), xc(3,k1), phia_new(k1)
        enddo
    !endif
    close(120) 
  
    call adhesion_tension_alternative
    !*******************************************************************!
    !                        CALCULATE NEW FIELD                        !
    !*******************************************************************!
    do k1 = 1, numnp
        distance   = xc(1,k1)
        wa_new(k1) = kapa * (phia_new(k1)- 1.d0) + Ufield(k1)
    enddo
    !*******************************************************************!
    !                     APPLY CONVERGENCE CRITERION                   !
    !*******************************************************************!
    error = 0.d00
    do k1 = 1, numnp
        error = max(error, dabs(wa_new(k1) - wa(k1)))
    enddo
    !*******************************************************************!
    !                     APPLY ANDERSON MIXING RULE                    !
    !*******************************************************************!
    do k1 = 1, numnp 
        wa(k1) = (1.d0-fraction) * wa(k1) + fraction * wa_new(k1)
    enddo

    if (mod(kk,1000).eq.0d0) then
        open(unit=123, file = 'field_compare.txt')
        do k1 = 1, numnp
            write(123,'(I10,2X,E16.9,2X,E16.9)') k1, wa(k1), wa_new(k1)
        enddo
    endif

    if (mod(kk,100).eq.0.d0) then
        write(6,*)  kk, adh_ten_alt, error
        write(iow,'(I10,2X,E16.9,2X,E16.9)') kk, adh_ten_alt, error
    endif

enddo!kk

if (pr_on==1) then 
    call qprint
endif

if (error.lt.max_error) then
    write(iow,'(/''Convergence of max error'',F16.9)') error
else
    write(iow,'(/''Convergence of '',I10, '' iterations'')') iterations
endif

write(iow,'(''-----------------------------------'')')
write(iow,'(''Adhesion tension (mN/m) '',E16.9)') adh_ten_alt
write(iow,'(''Partition function Q =  '',E16.9)') part_func
write(iow,'(''            n/n_bulk =  '',E16.9)') nch_per_area * dfloat(chainlen) / (rho_0 * lx  * 1.d-10)
!------------------------------------------------------------------------------------------------------------------!
end program FEM_3D
      