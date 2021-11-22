 subroutine mesh_io_3d
    use mdata
    use xdata
    use kcw
    implicit none
    
    
    character(len=15)::dummy
    integer ::idummy
    integer :: i,sdim,j!numnp,
    integer :: nmeltypes
    integer:: vtxnum,vtxel
    integer, allocatable::vtxelement(:,:)
    integer::vtxparnum, vtxparel
    integer, allocatable::vtxpar(:,:)
    integer,allocatable::vtxentity(:)
    integer:: edgnum,edgel
    integer, allocatable::edgelement(:,:)
    integer:: edgparnum,edgparel
    real,allocatable::edgpar(:,:)
    integer,allocatable::edgentity(:)
    !integer:: fcnum,fcel
    !integer, allocatable::fcelement(:,:)
    integer:: fcparnum,fcparel
    real,allocatable::fcpar(:,:)
    !integer,allocatable::fcentity(:)
    real,allocatable::temp3(:)
    real(8)::sum_f,Q

    
    
         open (unit=12, file='Fem_3D.mphtxt')
         open(unit=13, file='com.txt')
         do i=1,17
              read (12,'(a60)') dummy
         end do
         read (12,*) sdim
         read (12,*) numnp
     
         write(13,'(i9)') sdim
         write(13,'(i9)') numnp
         allocate (xc(sdim,numnp))
         do i=1,3
             read (12,'(a60)') dummy
         end do
! edw mporei na mpei ena format me to dimension                            
                             
         do i=1,numnp
              read(12,*)  (xc(j,i),j=1,sdim)
              write(13,'(f19.15,1x,f19.15,1x,f19.15)')(xc(j,i),j=1,sdim)
         end do
           
           
         lx=xc(1,numnp)
           
         open(55, file = 'mesh.txt')
         do i=1,numnp
           
        write(55,'(f21.15,1x,f21.15,1x,f21.15)')(xc(j,i),j=1,sdim) 
        end do
        close (55)
        read (12,'(a60)') dummy
!element types                 
        read (12,*)nmeltypes
        write (13,*) nmeltypes
                 
        do i=1,6
              read (12,'(a60)') dummy
        end do
                

                 
!vertex                 
!                 read (12,*) vtxnum
!                   read (12,*) vtxel
!                   allocate(vtxelement(vtxnum,vtxel))
!                   read (12,'(a60)') dummy
!                   do i=1, vtxel
!                       read(12,*) (vtxelement(j,i),j=1,vtxnum)
!                       write(13,*)(vtxelement(j,i),j=1,vtxnum)
!                   enddo
                 
                 
         call allocate_matrix(vtxnum,vtxel)
                 
         allocate (vtxelement(vtxnum,vtxel))
         call read_integer (vtxnum,vtxel,vtxelement)
                   
         read (12,'(a60)') dummy
                   
         call allocate_matrix(vtxparnum,vtxparel)
                 
         allocate(vtxpar(vtxparnum,vtxparel))
         call read_integer(vtxparnum,vtxparel,vtxpar)
!DEBUG
!                  read (12,'(a60)') dummy
!               if (vtxparel.ne.0)then
!                 do i=1,vtxparnum
!                     read(12,*) (vtxpar(j,i),j=1, vtxparel)
!                     write(13,*)(vtxpar(j,i),j=1, vtxparel)
!                  enddo
!                endif
!DEBUG
         read (12,'(a60)') dummy
         allocate(vtxentity(vtxel))
         call entity(vtxel,vtxentity)
         do i=1,8
              read (12,'(a60)') dummy
         end do
                    
!edge          
         call allocate_matrix(edgnum,edgel)
         allocate (edgelement(edgnum,edgel))
         call read_integer (edgnum,edgel,edgelement)
                 
                 
         read (12,'(a60)') dummy
         call allocate_matrix(edgparnum,edgparel)
         allocate(edgpar(edgparnum,edgparel))
         call read_real(edgparnum,edgparel,edgpar)
         read (12,'(a60)') dummy
         allocate(edgentity(edgel))
         call entity(edgel,edgentity)
         do i=1,9
              read (12,'(a60)') dummy
         end do
                    
               
                
                
!faces        
         call allocate_matrix(fcnum,fcel)
         allocate (fcelement(fcnum,fcel))
                    
         call read_integer (fcnum,fcel,fcelement)  
         read (12,'(a60)') dummy
! changes the values of face element so as to start from 1 instead from 0        
         fcelement=fcelement +1
         call allocate_matrix(fcparnum,fcparel)
         allocate (fcpar(fcparnum,fcparel))
         call read_real(fcparnum,fcparel,fcpar)
         read (12,'(a60)') dummy
         allocate(fcentity(fcel))
         call entity(fcel,fcentity)
! changes the values of  element entinty so as to start from 1 instead from 0        
         fcentity=fcentity +1         
         
         

                    
 !up/down pairs
         read (12,'(a60)') dummy
         read (12,* ) idummy
         do i=1,idummy
              read (12,'(a60)') dummy
         end do
                    
         do i=1,6
              read (12,'(a60)') dummy
         end do
!domain
                 
                   
         call allocate_matrix(dmnum,dmel)            
         allocate (ix(dmnum,dmel))
                    
         call read_integer (dmnum,dmel,ix)  
                    
         do i=1,dmel
              do j=1,dmnum
                  ix(j,i)=ix(j,i)+1
             end do
         end do 
         ndm =sdim
         nel=dmnum
         numel=dmel    
         allocate(temp3(numel))
         do i=1,numel
              temp3(i)=ix(7,i)
              ix(7,i)=ix(6,i)
              ix(6,i)=temp3(i)
         end do            
                    
                    
                    
                    
                    
         open (unit=77, file='inter.txt')
                    
         do i=1,dmel
              write(77,'(i9,10(2x,i9))')(ix(j,i),j=1,dmnum) 
         end do
         close (77)
                 
         
        
     !  open (file='out.txt',unit=iow)
       
        
!         do i=1,numnp
         
!         write(iow, '(e16.9)')xc(1,i)
!         end do 
     !  u=5.
         
         
                ! call spat_3d(u,sum_F,Q)
             !  call interpolation_3d(u,sum_F,Q)
                 !allocate(c_m%value(numnp))
                 !
                 !c_m%value(1)=7.98
                 !print*,'type success',c_m%value(1)
         allocate(c(numnp,numnp))
         allocate(k(numnp,numnp))
         allocate(w(numnp,numnp))
         allocate(rh(numnp,numnp))
         allocate(g(numnp,numnp))
                    
                
    return
end subroutine 
    
subroutine allocate_matrix(num,el)
    
    implicit none
    integer::num,el
    
         read(12,*)num
         read(12,*)el
    
end subroutine 
    
subroutine read_integer(num,el,element)
    
    implicit none
    integer:: num,el,i,j
    character (len=60)::dum
    integer::element(num,el)

         read (12,'(a60)') dum
         if(el.ne.0)then
              if(num.eq.1)then
                   do i=1,el
                        read(12,*)element(1,i)
                        write(13,*)element(1,i)
                   end do
              else 
                  do i=1,el
                        read(12,*)(element(j,i),j=1,num)
                        write(13,*)(element(j,i),j=1,num)
                 end do
              end if
         end if
end subroutine
    
subroutine entity(el,endity)
    integer::num,el,equal,i
    integer::endity(el)
    character (len=60)::dum
    read(12,*)equal
    print*,equal,el
    read (12,'(a60)') dum
             if (equal.eq.el)then
                 do i=1,el
                     read(12,*) endity(i)
                     write(13,*)endity(i)
                 end do
             else 
                 print*, "error in entity"
             end if
end subroutine
subroutine read_real(num,el,element)
    
    implicit none
    integer:: num,el,i,j
    character (len=60)::dum
    real::element(num,el)
   
   
     read (12,'(a60)') dum
     if(el.ne.0)then
     
              if(num.eq.1)then
                    do i=1,el
                         read(12,*)element(1,i)
                         write(13,*)element(1,i)
                 end do
         else 
             do i=1,el
                  read(12,*)(element(j,i),j=1,num)
                write(13,*)(element(j,i),j=1,num)
             end do
        end if
     end if
    end subroutine
  