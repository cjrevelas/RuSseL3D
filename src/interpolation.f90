
subroutine interpolation_3d(u_spat,sum_out,Q)

    use xdata
    use mdata
    
    implicit none

    real*8 Q
    real*8 ul(nel),xl(ndm,nel),xl1(ndm,nel)
    integer lint
    real*8 shp(4,11), sv(5,11), xsj,xsj1,u_spat(numnp)
    real*8 uqp, sumel,sumel1,sum_out,sum_out1, volel,volel1, vol,vol1
    real*8 time
    integer i, j, ii, l,n,ne4,k
      ne4=4

!Calculation of integral in space
!==================================================================!
!Initialize accumulator for integral of the solution
         sum_out = 0.d00
!Initialize accumulator for total domain volume
         vol = 0.d00
!Loop over all elements
!------------------------------------------------------------------!
         do n = 1, numel
!  Loop over all nodes of current element
               do i = 1, ne4
!          Find global index of node
                   ii = ix(i,n)                       
                   do j = 1, ndm
                        xl(j,i) = xc(j,ii)
                   end do!j       
!      Copy value of solution at current node from global to local array
                   ul(i) = u_spat(ii) 
               end do!i
!      Set up for  quadrature. 
            
              !lint =1 
             l=3
             ! call gauss_1d(lint,sv)
             call gauss_3d(l,lint,sv)
              sumel = 0.d00

              volel = 0.d00

!      Loop over all quadrature points in element
              do l = 1,lint
                  
                  ! call   shp1d(sv(1,l),xl,shp,ndm,nel,xsj)
                call   tetshp( sv(1,l), xl, ndm, ne4, xsj, shp )
             !  call     tetshp10( sv(1,l), xl, ndm, xsj, shp)
                  ! xsj = xsj*sv(2,l)
                   xsj = xsj*sv(5,l)
                   uqp   = 0.0d0
                   
                   do j = 1, ne4
                        uqp = uqp + shp(4,j)*ul(j)
                   end do!j
                   
                   sumel = sumel + uqp * xsj
                   volel = volel + xsj
                  
              end do!l
              sum_out = sum_out + sumel
              vol = vol + volel
              


         end do!n  End loop over elements
         !stop
         Q = sum_out/vol

!DEBUG
!------
  write(iow,'('' SUM '',E16.9,'' volume '',E16.9,'' mean value  of solution '',E16.9)') sum_out, vol, q
!  write(99,'(E16.9,2X,E16.9)') ttim, sum/vol
!DEBUG
    return
end


      
