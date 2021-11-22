 subroutine Matrix_assemble(C_M,K_M,W_M)
    
     use mdata
     use xdata
     
    implicit none
    
    integer :: lint
    integer :: m_1,n_1
    real*8 ul(ndm,nel,7),xl(ndm,nel)
    real*8 shp(4,11), sv(5,11), xsj
    integer i, j,ii,l,n,m,nn
    double precision :: C_M(numnp,numnp)
    double precision :: K_M(numnp,numnp)
    double precision :: W_M(numnp,numnp)
    integer,dimension (nel) :: gln
          
            
              
         K_M=0.
         C_M=0.
         W_M=0.
         do nn = 1, numel
!          1.Loop over all nodes of current element      
!          2. Find global index of node     
!          3.Copy coordinates from global array to local array
!          concerning current element      
              do i = 1, nel
                  gln(i)= ix(i,nn)
                  do j = 1, ndm
                        xl(j,i) = xc(j, gln(i))
                   end do!j
              end do!i
!      Set up for gauss quadrature. 
             !lint =1 
              l=3
            ! call gauss_1d(lint,sv)
            
            call gauss_3d(l,lint,sv)
!      Loop over all points used for quadrature in the considered element.
!      nint:  Total number of quadrature points in each element 
!             (11 for second-order tetrahedra in thermal problem).

!      Loop over all quadrature points in element
!          sv(2,l)------->weight
!          shp(2,nel): Shape functions and derivatives at s;
!          shp(1,1 to nel): derivatives of shape functions
!          shp(2,1 to nel): shape functions             
              do l = 1,lint
                  ! call shp1d(sv(1,l),xl,shp,ndm,nel,xsj)
                  call tetshp10( sv(1,l), xl, ndm, xsj, shp)
    
                  do m = 1,nel
                       m_1 = gln(m)
                       do n = 1,nel
                             n_1 = gln(n)
                            ! C_M(m_1,n_1)=C_M(m_1,n_1) + shp(2,n)*shp(2,m)*xsj*sv(2,l)
                            C_M(m_1,n_1)=C_M(m_1,n_1)+shp(4,n)*shp(4,m)*xsj*sv(5,l)
                            ! K_M(m_1,n_1)=K_M(m_1,n_1)+(Rgyr**2.)*shp(1,n)*shp(1,m)*xsj*sv(2,l)
                             K_M(m_1,n_1)=K_M(m_1,n_1)+&
										(Rgyr**2.)*(shp(1,n)*shp(1,m)+shp(2,n)*shp(2,m)+shp(3,n)*shp(3,m))*xsj*sv(5,l)
                             W_M(m_1,n_1)=W_M(m_1,n_1) + wa(n_1)*shp(4,n)*shp(4,m)*xsj*sv(5,l)
                        enddo !n
                   enddo !m
              end do!l
         end do

!DEBUG             
             !open (file ='k.txt',unit=23)
             !open(file='c.txt',unit=24)
             !open(file='w.txt',unit=25)
             !
             !do i=1,numnp
             !   write (23,'(<numnp>e19.9)') (k_m(i,n_1),n_1=1,numnp)
             !   write (24,'(<numnp>e19.9)') (c_m(i,n_1),n_1=1,numnp)
             !    write (25,'(<numnp>e19.9)') (w_m(i,n_1),n_1=1,numnp)
             !end do 
!DEBUG           
    return
end subroutine

    
