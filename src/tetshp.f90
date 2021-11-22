
      subroutine tetshp( xi, xl, ndm, nel, xsj, shp )

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2014: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/11/2006
!       1. Limit number of element nodes for order = 1      11/08/2007
!       2. Reprogram using xi(1) = 1.- xi(2)- xi(3) - xi(4) 12/08/2007
!          Change jacobian xsj to volume of tetrahedron
!       3. Add 14- & 15-node shape functions with face dof  31/08/2007
!       4. Add heirarchic bubble for 4 node element (-1)    21/09/2007
!       5. Add heirarchic bubble for 10/14 node element(-2) 26/09/2007
!       6. Change argument from 'order' to 'nel'            05/11/2007
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Compute 3-d tetrahedral element shape
!               functions and their derivatives w/r x,y,z

!      Inputs:
!         xi(4)     - Natural volume coordinates of point
!         xl(ndm,*) - Nodal coordinates for element
!         ndm       - Spatial dimension of mesh
!         nel       - Number of nodes/ element
!                   -  4-nodes     = linear
!                   - 10-/11-nodes = quadratic + bubble
!                   - 14-/15-nodes = quadratic + face + bubble

!      Outputs:
!         xsj       - Jacobian determinant at point
!         shp(4,*)  - Shape functions and derivatives at point
!                     shp(1,i) = dN_i/dx
!                     shp(2,i) = dN_i/dy
!                     shp(3,i) = dN_i/dz
!                     shp(4,i) =  N_i
!-----[--.----+----.----+----.-----------------------------------------]
      implicit none

!      include 'iofile.h'
!      include 'pconstant.h'

      integer  ndm, nel, nnel, i, j, k
      real*8   xsj, detr, xi(4), xl(ndm,*), shp(4,*), xs(3,3),xsi(3,3)
      real*8   a(3,3), te(3)

      save

!     Linear shape functions and their derivatives

      if(abs(nel).eq.4) then
        shp(1, 1) =-1.d0
        shp(1, 2) = 0.d0
        shp(1, 3) = 0.d0
        shp(1, 4) = 1.d0

        shp(2, 1) =-1.d0
        shp(2, 2) = 1.d0
        shp(2, 3) = 0.d0
        shp(2, 4) = 0.d0

        shp(3, 1) =-1.d0
        shp(3, 2) = 0.d0
        shp(3, 3) = 1.d0
        shp(3, 4) = 0.d0

        shp(4, 1) = xi(1)
        shp(4, 2) = xi(2)
        shp(4, 3) = xi(3)
        shp(4, 4) = xi(4)
        nnel      = 4

!     Quadratic shape functions and derivatives

      elseif(abs(nel).eq.10 .or. abs(nel).eq.11  .or.&
            abs(nel).eq.14 .or. abs(nel).eq.15) then
          
      shp(1, 1) = 4.d0*xi(1) - 1.d0
      shp(1, 2) = 0.d0
      shp(1, 3) = 0.d0
      shp(1, 4) =-4.d0*xi(4) + 1.d0
      shp(1, 5) = 4.d0*xi(2)
      shp(1, 6) = 0.d0
      shp(1, 7) = 4.d0*xi(3)
      shp(1, 8) = 4.d0*(xi(4) - xi(1))
      shp(1, 9) =-4.d0*xi(2)
      shp(1,10) =-4.d0*xi(3)

      shp(2, 1) = 0.d0
      shp(2, 2) = 4.d0*xi(2) - 1.d0
      shp(2, 3) = 0.d0
      shp(2, 4) =-4.d0*xi(4) + 1.d0
      shp(2, 5) = 4.d0*xi(1)
      shp(2, 6) = 4.d0*xi(3)
      shp(2, 7) = 0.d0
      shp(2, 8) =-4.d0*xi(1)
      shp(2, 9) = 4.d0*(xi(4) - xi(2))
      shp(2,10) =-4.d0*xi(3)

      shp(3, 1) = 0.d0
      shp(3, 2) = 0.d0
      shp(3, 3) = 4.d0*xi(3) - 1.d0
      shp(3, 4) =-4.d0*xi(4) + 1.d0
      shp(3, 5) = 0.d0
      shp(3, 6) = 4.d0*xi(2)
      shp(3, 7) = 4.d0*xi(1)
      shp(3, 8) =-4.d0*xi(1)
      shp(3, 9) =-4.d0*xi(2)
      shp(3,10) = 4.d0*(xi(4) - xi(3))
!
!   Compute shape functions and store in shp(4,*)

      shp(4, 1) = (2.d0*xi(1) - 1.d0) * xi(1)
      shp(4, 2) = (2.d0*xi(2) - 1.d0) * xi(2)
      shp(4, 3) = (2.d0*xi(3) - 1.d0) * xi(3)
      shp(4, 4) = (2.d0*xi(4) - 1.d0) * xi(4)
      shp(4, 5) = 4.d0 * xi(1) * xi(2)
      shp(4, 6) = 4.d0 * xi(2) * xi(3)
      shp(4, 7) = 4.d0 * xi(3) * xi(1)
      shp(4, 8) = 4.d0 * xi(4) * xi(1)
      shp(4, 9) = 4.d0 * xi(4) * xi(2)
      shp(4,10) = 4.d0 * xi(4) * xi(3)
          
          
          
          
          
          
          
          
          
          
          
          
          
          

        !shp(1, 1) =-4.d0*xi(1) + 1.d0
        !shp(1, 2) = 4.d0*xi(2) - 1.d0
        !shp(1, 3) = 0.d0
        !shp(1, 4) = 0.d0
        !shp(1, 5) = 4.d0*(xi(1) - xi(2))
        !shp(1, 6) = 4.d0*xi(3)
        !shp(1, 7) =-4.d0*xi(3)
        !shp(1, 8) =-4.d0*xi(4)
        !shp(1, 9) = 4.d0*xi(4)
        !shp(1,10) = 0.d0
        !
        !shp(2, 1) =-4.d0*xi(1) + 1.d0
        !shp(2, 2) = 0.d0
        !shp(2, 3) = 4.d0*xi(3) - 1.d0
        !shp(2, 4) = 0.d0
        !shp(2, 5) =-4.d0*xi(2)
        !shp(2, 6) = 4.d0*xi(2)
        !shp(2, 7) = 4.d0*(xi(1) - xi(3))
        !shp(2, 8) =-4.d0*xi(4)
        !shp(2, 9) =-0.d0
        !shp(2,10) = 4.d0*xi(4)
        !
        !shp(3, 1) =-4.d0*xi(1) + 1.d0
        !shp(3, 2) = 0.d0
        !shp(3, 3) = 0.d0
        !shp(3, 4) = 4.d0*xi(4) - 1.d0
        !shp(3, 5) =-4.d0*xi(2)
        !shp(3, 6) = 0.d0
        !shp(3, 7) =-4.d0*xi(3)
        !shp(3, 8) = 4.d0*(xi(1) - xi(4))
        !shp(3, 9) = 4.d0*xi(2)
        !shp(3,10) = 4.d0*xi(3)
        !
        !shp(4, 1) = xi(1)*(2.d0*xi(1) - 1.d0)
        !shp(4, 2) = xi(2)*(2.d0*xi(2) - 1.d0)
        !shp(4, 3) = xi(3)*(2.d0*xi(3) - 1.d0)
        !shp(4, 4) = xi(4)*(2.d0*xi(4) - 1.d0)
        !shp(4, 5) = 4.d0*xi(1)*xi(2)
        !shp(4, 6) = 4.d0*xi(2)*xi(3)
        !shp(4, 7) = 4.d0*xi(3)*xi(1)
        !shp(4, 8) = 4.d0*xi(1)*xi(4)
        !shp(4, 9) = 4.d0*xi(2)*xi(4)
        !shp(4,10) = 4.d0*xi(3)*xi(4)

!       14- & 15-node tetrahedron: Set face and internal node

        if(nel.eq.14 .or. nel.eq.15) then

!         Mid-face shape functions

          shp(1,11) =  27.0d0*(xi(1) - xi(2))*xi(4)
          shp(1,12) =  27.0d0*xi(3)*xi(4)
          shp(1,13) = -27.0d0*xi(3)*xi(4)
          shp(1,14) =  27.0d0*(xi(1) - xi(2))*xi(3)

          shp(2,11) = -27.0d0*xi(2)*xi(4)
          shp(2,12) =  27.0d0*xi(2)*xi(4)
          shp(2,13) =  27.0d0*(xi(1) - xi(3))*xi(4)
          shp(2,14) =  27.0d0*(xi(1) - xi(3))*xi(2)

          shp(3,11) =  27.0d0*(xi(1) - xi(4))*xi(2)
          shp(3,12) =  27.0d0*xi(2)*xi(3)
          shp(3,13) =  27.0d0*(xi(1) - xi(4))*xi(3)
          shp(3,14) = -27.0d0*xi(2)*xi(3)

          shp(4,11) =  27.0d0*xi(1)*xi(2)*xi(4)
          shp(4,12) =  27.0d0*xi(2)*xi(3)*xi(4)
          shp(4,13) =  27.0d0*xi(3)*xi(1)*xi(4)
          shp(4,14) =  27.0d0*xi(1)*xi(2)*xi(3)

!         Correct vertex and mid-edge shape functions

          do j = 1,4
            shp(j, 1) = shp(j, 1)&
                     +( 1./9.) *(shp(j,13) + shp(j,11) + shp(j,14))
            shp(j, 2) = shp(j, 2)&
                     + ( 1./9.) *(shp(j,11) + shp(j,12) + shp(j,14))
            shp(j, 3) = shp(j, 3)& 
                    +( 1./9.) *(shp(j,12) + shp(j,13) + shp(j,14))
            shp(j, 4) = shp(j, 4)&
                     + ( 1./9.) *(shp(j,11) + shp(j,12) + shp(j,13))
            shp(j, 5) = shp(j, 5) - (4./9.)*(shp(j,11) + shp(j,14))
            shp(j, 6) = shp(j, 6) - (4./9.)*(shp(j,12) + shp(j,14))
            shp(j, 7) = shp(j, 7) - (4./9.)*(shp(j,13) + shp(j,14))
            shp(j, 8) = shp(j, 8) - (4./9.)*(shp(j,13) + shp(j,11))
            shp(j, 9) = shp(j, 9) - (4./9.)*(shp(j,11) + shp(j,12))
            shp(j,10) = shp(j,10) - (4./9.)*(shp(j,12) + shp(j,13))
          end do ! j
        endif

!       Interior node functions for 11 and 15-node element

        if(nel.eq.11 .or. nel.eq.15) then

          shp(1,nel) = 256.d0*(xi(1) - xi(2))*xi(3)*xi(4)
          shp(2,nel) = 256.d0*(xi(1) - xi(3))*xi(2)*xi(4)
          shp(3,nel) = 256.d0*(xi(1) - xi(4))*xi(2)*xi(3)
          shp(4,nel) = 256.d0*xi(1)*xi(2)*xi(3)*xi(4)

!         Correct vertex shape functions for interior values

          do j = 1,4
            do i = 1,4
              shp(i,j) = shp(i,j) + 0.125d0*shp(i,nel)
            end do ! i
          end do ! j

!         Correct mid-edge shape functions for interior values

          do j = 5,10
            do i = 1,4
              shp(i,j) = shp(i,j) - 0.250d0*shp(i,nel)
            end do ! i
          end do ! j
        endif

        nnel = nel

!     Error - Higher than quadratic not coded

      else

!        write(iow,2000) nel
        write(  *,2000) nel
!        call plstop()

      endif

!     Compute jacobian matrix

      do i = 1,3
        do j = 1,3
          xs(j,i) = 0.0d0
          do k = 1,nnel
            xs(j,i) = xs(j,i) + xl(j,k)*shp(i,k)
          end do ! k
        end do ! i
      end do ! j

!     Compute inverse of jacobian matrix

      xsi(1,1) = xs(2,2)*xs(3,3) - xs(3,2)*xs(2,3)
      xsi(1,2) = xs(3,2)*xs(1,3) - xs(1,2)*xs(3,3)
      xsi(1,3) = xs(1,2)*xs(2,3) - xs(2,2)*xs(1,3)

      xsi(2,1) = xs(2,3)*xs(3,1) - xs(3,3)*xs(2,1)
      xsi(2,2) = xs(3,3)*xs(1,1) - xs(1,3)*xs(3,1)
      xsi(2,3) = xs(1,3)*xs(2,1) - xs(2,3)*xs(1,1)

      xsi(3,1) = xs(2,1)*xs(3,2) - xs(3,1)*xs(2,2)
      xsi(3,2) = xs(3,1)*xs(1,2) - xs(1,1)*xs(3,2)
      xsi(3,3) = xs(1,1)*xs(2,2) - xs(2,1)*xs(1,2)

      do i = 1,3
        do j = 1,3
          a(i,j) = xs(i,1)*xsi(1,j)+xs(i,2)*xsi(2,j)+xs(i,3)*xsi(3,j)
        end do ! j
      end do ! i

!     Compute jacobian determinant

      xsj = xs(1,1)*xsi(1,1) + xs(1,2)*xsi(2,1) + xs(1,3)*xsi(3,1)

      if(xsj.ne.0.0d0) then
        detr = 1.d0/xsj
        !change by tolis
        xsj  = xsj!*(1./6.)
        !change by tolis
      else
!        write(iow,*) ' TETSHP: Determinant =',xsj
        detr = 1.d0
      endif

!     Compute jacobian inverse

      do j = 1,3
        do i = 1,3
          xs(i,j) = xsi(i,j)*detr
        end do ! i
      end do ! j

!     Heirarchic interior node function for 4-node

      if(nel .eq. -4) then
        shp(1,5) = 256.d0*(xi(1) - xi(2))*xi(3)*xi(4)
        shp(2,5) = 256.d0*(xi(1) - xi(3))*xi(2)*xi(4)
        shp(3,5) = 256.d0*(xi(1) - xi(4))*xi(2)*xi(3)
        shp(4,5) = 256.d0*xi(1)*xi(2)*xi(3)*xi(4)
        nnel     = 5

!     Hierarchical bubble function for 10 and 14-node element

      elseif(nel.eq.-10 .or. nel.eq.-14) then

        if(abs(nel).eq.10) then
          nnel = 11
        else
          nnel = 15
        endif
        shp(1,nnel) = 256.d0*(xi(1) - xi(2))*xi(3)*xi(4)
        shp(2,nnel) = 256.d0*(xi(1) - xi(3))*xi(2)*xi(4)
        shp(3,nnel) = 256.d0*(xi(1) - xi(4))*xi(2)*xi(3)
        shp(4,nnel) = 256.d0*xi(1)*xi(2)*xi(3)*xi(4)

      endif

!     Compute shape function derivatives

      do k = 1,nnel
        do i = 1,3
          te(i) = shp(1,k)*xs(1,i) + shp(2,k)*xs(2,i) + shp(3,k)*xs(3,i)
        end do ! i
        do i = 1,3
          shp(i,k) = te(i)
        end do ! i
      end do ! k

!     Format

2000  format(/' *ERROR* TETSHP not coded for nel =',i4)

      end
