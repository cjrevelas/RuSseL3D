subroutine edwards_free_film_fem
    
  
    use xdata
    use mdata 
    use kcw
    
    
    
    implicit none
     integer::i,j,f,idummy,fc,i9
     double precision, allocatable ::u1(:)


!===================================================================!
!                                                                   !
!            FREE Edwards Diffusion Equation Solution               !
!                                                                   !
!===================================================================!
!-------------------------------------------------------------------!
!                                                                   !
!........Initial onditions for homopolymer                          !
!........diffusion equations                                        !
!                                                                   !
!-------------------------------------------------------------------!
         
!-------------------------------------------------------------------!
!........The potential is set zero at x=infinity                    !
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
!...........In qf(*,1) is stored the propagator q(x,s)              !
!...........for the !urrent value of s, where s in!reases           !
!...........from 0 to 1 as we go from !hain start to !hain end.     !
!...........In qf(*,2) is formed the propagator q(x,s+ds)           !
!...........for the !urrent value of s.                             !
!...........Note that this !ode has resulted from adaptation        !
!...........of a previous !ode for !opolymers to homopolymers       !
!........... q2(*,1) in NOT USED                                    !
!                                                                   !
!...........Propagator values for all s, from !hain start           !
!...........to !hain end, are stored in qf_final(*,*)               !
!...........Propagator values for all s, from !hain end             !
!...........to !hain start, are stored in q2_final(*,*)             !
!...........For the homopolymer treated here, qf_final              !
!...........and q2_final are identical.                             !
!!..................................................................!
!........Initial Conditions                                         !
!.................................................................. !
!...........Initial value of propagator, q(x,0) = 1.0 for all x     !
!.............The initial values stored to qf_final for s=0         !
     
         do i1 = 1,numnp
              qf(i1,1) = 1.d0
              qf_final(i1,1) = 1.d0
         end do!i1
         
!-------------------------------------------------------------------!
!........Boundary Conditions                                        !
!                                                                   !
!........Absorbing boundary Condition at x = 0                      !
!........Bulk Conditions prevailing at   x = lx                     !
!-------------------------------------------------------------------!
!             qf(1,1)    = 0.d0
!-------------------------------------------------------------------!
!........Implement Semi implicit method                             !
!........in order to solve the diffusion equation                   !
!........for  qf                                                    !
!........In order to save memory propagator                         !
!........will be Calculated locally for one timestep                !
!........end the volume fra!tions will be calculated                !
!........on the fly                                                 !
!-------------------------------------------------------------------!
!........"Time" here corresponds to moving along the contour        !
!-------------------------------------------------------------------!
allocate(u1(numnp))
         do time_step = 2, ns+1
              
              
              
!------------------------------------------------------------------ !
!                   Finite element Method                           !
!-------------------------------------------------------------------!

              g=0.
              rh=0.
              do i =1,numnp
                   do j =1,numnp
!                        g(i,j)  =   c(i,j) +(ds/2.d0)* (k(i,j)+w(i,j))
 !                       rh(i,j) = c(i,j) -(ds/2.d0)* (k(i,j)+w(i,j))
                         g(i,j)  =   c(i,j) +(ds)* (k(i,j)+w(i,j))
                         rh(i,j) = c(i,j) !-(ds/2.d0)* (k(i,j)+w(i,j))
                   enddo !j
              enddo !i
             
             
!              do i = 1,numel+1
!                   bdiag(i)=g(i,i)
!              end do!i
!              do i = 2,numel+1
!                   adiag(i)=g(i,i-1)
!              end do!i
!              do i = 1,numel
!                   cdiag(i)=g(i,i+1)
!              end do!i
!OLD trigad Creation            
              !do i = 1,numnp
              !  do j = 1,numnp
              !
              !      if ( j.eq.i )   bdiag(j) = g(i,j)
              !      if ((j.ne.1).and.(j.eq.(i)) )  adiag(j)= g(i,j)          
              !      if ((i.lt.numnp).and.(j.eq.(i+1))) cdiag(j) = g(i,j)
              !
              !  enddo !j
              !enddo !i
!OLD trigad Creation              

              rdiag=0.
              do i =1,numnp
                   do j =1,numnp
                        rdiag(i) = rdiag(i) + rh(i,j)*qf(j,1)
                   enddo !j
              enddo !i
!AYTH H LOOPA EINAI AHDIA
              
              
         
        do j=1,fcel
            if ((fcentity(j)==3).or.(fcentity(j)==4) )then
                do i=1,fcnum
                    idummy= fcelement(i,j)
                    do f=1,numnp
                        g(idummy ,f)=0
                    end do 
                    g(idummy,idummy)=1.
                    rdiag(idummy)=0.
                end do 
            end if 
        end do
            


              !bdiag(1) = 1.d0
              !cdiag(1) = 0.d0
              !rdiag(1) = 0.d0
              

!         call tridag(adiag,bdiag,cdiag,rdiag,u,numnp)
            non_zero=0
            
            
        if (time_step==2)then
        i9=0
              do i=1,numnp
                  do j=1,numnp
                      if(abs(g(i,j))>1.e-8) then
                          non_zero=non_zero+1
                      end if 
                  end do 
              end do 
              allocate (c_m%value(non_zero))
              allocate (c_m%col(non_zero))
              allocate (c_m%row(non_zero))
        end if 
       
              
                non_zero=0
              do i=1,numnp
                  do j=1,numnp
                      if(abs(g(i,j))>1.e-8) then
                          non_zero=non_zero+1
                          c_m%value(non_zero)=g(i,j)
                          c_m%row(non_zero)=i
                          c_m%col(non_zero)=j
                      end if 
                  end do
              end do 
              
                 
             open (file ='cm.txt',unit=99)
                do i=1,non_zero
                   write (99,'(i6,i6,e16.9)' ) c_m%row(i),c_m%col(i) ,c_m%value(i)   
                end do 
             close(99)
             
                call  MUMPS_SUB(numnp,rdiag,u1,i9)
                
!DEBUG              
                if(time_step==2)then
                i9=45
                    open(file='x1.txt',unit=23)
                   do i1=1,numnp
                         write(23,'(i4,2x,f10.4)')i1,rdiag(i1)
                    end do 
                end if
!DEBUG              
                
                
                do i1 = 1,numnp
                  qf(i1,2) = u1(i1)
               end do!i1
!DEBUG             
              
!            open (file ='g.txt',unit=93)
!             open(file='rh.txt',unit=94)
!            open(file='r.txt',unit=95)
             
!             do i=1,numnp
!                  write (93,'(<numnp>e19.9)') (g(i,j),j=1,numnp)
!                  write (94,'(<numnp>e19.9)') (rh(i,j),j=1,numnp)
!                  write (95,'(e19.9)') rdiag(i)
!             end do
                 
!DEBUG                 
!-------------------------------------------------------------------!
!................Save the propagators                               !
!................in  qf_final(x,s) [s from start to end]            !
!-------------------------------------------------------------------!
!................Copy q(*,s+ds) into q(*,s) to prepare for          !
!................next step in s.                                    !
!-------------------------------------------------------------------!    
              
              do i1 = 1,numnp
                   qf_final(i1,time_step) = qf(i1,2)
                   qf(i1,1) = qf(i1,2)
              end do !i1

         end do
!-------------------------------------------------------------------!         
!---Solution of Edwards Difusion equation for this iteration--------!
!-------------------------------------------------------------------!          
!DEBUG
        !open (file='qall.txt',unit=96)
        ! write(96,'(<ns+1>i19)')( i1,i1=0,ns)
        !do i1=0,nx
        !   write(96,'(<ns+1>e19.9)') (qf_final(i1,time_step),time_step=0,ns)
        ! end do 


!DEBUG                  
            
              
              
    return 
end
