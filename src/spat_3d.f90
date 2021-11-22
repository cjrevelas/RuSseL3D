subroutine spat_3d(u_spat,sum_out,Q)

    use xdata
    use mdata
    
    implicit none

    real*8 ttim, Q
    real*8 ul(nel),xl(ndm,nel)
    integer lint
    real*8 shp(4,11), sv(5,11), xsj,u_spat(numnp)
    real*8 uqp, sumel,sum_out, volel, vol
    real*8 time
    integer i, j, ii, l,n
      

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
               do i = 1, nel
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
                  
                  !DEBUG
                    call tetshp10( sv(1,l), xl, ndm, xsj, shp)
                     !call tetshp( sv(1,l), xl, ndm,nel, xsj, shp)
                     
                  !DEBUG
                  ! xsj = xsj*sv(2,l)
                   xsj = xsj*sv(5,l)
                   uqp   = 0.0d0
                   
                   do j = 1, nel
                        uqp = uqp + shp(4,j)*ul(j)
                   end do!j
                   
                   sumel = sumel + uqp * xsj
                   volel = volel + xsj
                  
              end do!l
              sum_out = sum_out + sumel
              vol = vol + volel

         end do!n  End loop over elements
         Q = sum_out/vol

!DEBUG
!------
  write(iow,'('' SUM '',E16.9,'' volume '',E16.9,'' mean value  of solution '',E16.9)') sum_out, vol, q
!  write(99,'(E16.9,2X,E16.9)') ttim, sum/vol
!DEBUG
    return
end


      