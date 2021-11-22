
!*** Begin change by DNT
      subroutine tetshp10( xi, xl, ndm, xsj, shp )

!-----[--.----+----.----+----.-----------------------------------------]
!    Purpose: Compute 3-d 10-node (second order) tetrahedral element 
!             shape functions and their derivatives w/r r_1, r_2, r_3.

!    Inputs:
!       xi(4)     - Natural volume coordinates of point
!       xl(ndm,*) - Nodal coordinates for element
!       ndm       - Spatial dimension of mesh

!    Outputs:

!       shp(4,*)  - Shape functions and derivatives with respect
!                   to nodal coordinates at point
!                   shp(1,k) =  partial N^k/ partial r_1
!                   shp(2,k) =  partial N^k/ partial r_2
!                   shp(3,k) =  partial N^k/ partial r_3
!                   shp(4,k) =  N^kc
!       xs(3,3)   - transpose of Jacobian matrix 
!                                                       
!                   xs(i,j) = partial r_j/partial xi_i =

!                             10
!                             Sum xl(j,k) (partial N^k/partial xi_i)
!                             k=1 

!                   with (r_1, r_2, r_3) being the coordinates of the
!                   considered point in space,

!                   (xi_1, xi_2, xi_3) being the natural coordinates
!                   of the same point in the element where it belongs

!                   and xl(j,k) being the j-th coordinate of node k
!                   connected to the element in the laboratory frame
!                   [k=1(1)10, j=1(1)3]

!                              10
!                   Note r_j = Sum xl(j,k) N^k(xi_1, xi_2, xi_3)
!                              k=1

!                   and, for any function f(r),
!                                           ~
!                          10
!                   f(r) = Sum f_k N^k(xi_1, xi_2, xi_3)
!                     ~    k=1  

!                   with f_k being the values of the function at the 
!                   nodes of the element where r belongs
!                                              ~

!       xsj       - Jacobian determinant at point, det(xs)
!Auxiliary

!       shpxi(3,10)- Derivatives of shape functions with respect
!                   to natural coordinates at point
!                   shpxi(1,k) =  partial N^k/ partial xi_1
!                   shpxi(2,k) =  partial N^k/ partial xi_2
!                   shpxi(3,k) =  partial N^k/ partial xi_3c
!       xsinv(3,3)- Jacobian inverse
!-----[--.----+----.----+----.-----------------------------------------]
!       use  iofile
!        use  pointer1
!       use  comblk
        implicit  none

      
      

      integer   i, j, k
      integer   ndm, nel

      real*8    xsj
      real*8    xi(4), xl(ndm,*), shp(4,*)
      real*8    xs(3,3)      

      real*8    shpxi(3,10), xsinv(3,3)

!   real*8    prod(3,3)
!
      save
!  
!   Compute shape function derivatives with respect to natural coordinates

!   shpxi(i,k) = partial N^k/partial xi_i            (i=1,2,3)

!   with N^k being the shape function for node k and
!   xi_i being the ith component of the natural coordinate vector.
!   Commands are the same as used in routine tjac3d to form array shp.

 
      xi(4) = 1.d0 - xi(1) - xi(2) - xi(3)

      shpxi(1, 1) = 4.d0*xi(1) - 1.d0
      shpxi(1, 2) = 0.d0
      shpxi(1, 3) = 0.d0
      shpxi(1, 4) =-4.d0*xi(4) + 1.d0
      shpxi(1, 5) = 4.d0*xi(2)
      shpxi(1, 6) = 0.d0
      shpxi(1, 7) = 4.d0*xi(3)
      shpxi(1, 8) = 4.d0*(xi(4) - xi(1))
      shpxi(1, 9) =-4.d0*xi(2)
      shpxi(1,10) =-4.d0*xi(3)

      shpxi(2, 1) = 0.d0
      shpxi(2, 2) = 4.d0*xi(2) - 1.d0
      shpxi(2, 3) = 0.d0
      shpxi(2, 4) =-4.d0*xi(4) + 1.d0
      shpxi(2, 5) = 4.d0*xi(1)
      shpxi(2, 6) = 4.d0*xi(3)
      shpxi(2, 7) = 0.d0
      shpxi(2, 8) =-4.d0*xi(1)
      shpxi(2, 9) = 4.d0*(xi(4) - xi(2))
      shpxi(2,10) =-4.d0*xi(3)

      shpxi(3, 1) = 0.d0
      shpxi(3, 2) = 0.d0
      shpxi(3, 3) = 4.d0*xi(3) - 1.d0
      shpxi(3, 4) =-4.d0*xi(4) + 1.d0
      shpxi(3, 5) = 0.d0
      shpxi(3, 6) = 4.d0*xi(2)
      shpxi(3, 7) = 4.d0*xi(1)
      shpxi(3, 8) =-4.d0*xi(1)
      shpxi(3, 9) =-4.d0*xi(2)
      shpxi(3,10) = 4.d0*(xi(4) - xi(3))
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
!
!   Compute transpose of Jacobian matrix
!
!                       partial r_j    10          partial N^k
!     xs(i,j) = J_ji = ------------- = Sum xl(j,k) ------------
!                       partial xi_i   k=1         partial xi_i
!
!   i,j in {1,2,3}
!   Minus sign is due to the fact that the 
!   laboratory system and the local system of 
!   natural coordinates of the element
!  have opposite handedness.

      do i = 1,3
        do j = 1,3
          xs(i,j) = 0.0d0
          do k = 1,10
            xs(i,j) = xs(i,j) - shpxi(i,k)*xl(j,k)
          end do
        end do
      end do

!  Compute jacobian determinant

      xsj  = xs(1,1)*(xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2))&
          + xs(1,2)*(xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3))&
          + xs(1,3)*(xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1))

!  Compute Jacobian inverse: 

!          xsinv(i,j) = (J^-1)_ij = partial xi_i/partial r_j

      xsinv(1,1) = xs(2,2)*xs(3,3) - xs(3,2)*xs(2,3)
      xsinv(2,1) = xs(3,2)*xs(1,3) - xs(1,2)*xs(3,3)
      xsinv(3,1) = xs(1,2)*xs(2,3) - xs(2,2)*xs(1,3)

      xsinv(1,2) = xs(3,1)*xs(2,3) - xs(2,1)*xs(3,3)
      xsinv(2,2) = xs(1,1)*xs(3,3) - xs(3,1)*xs(1,3)
      xsinv(3,2) = xs(2,1)*xs(1,3) - xs(1,1)*xs(2,3)

      xsinv(1,3) = xs(2,1)*xs(3,2) - xs(3,1)*xs(2,2)
      xsinv(2,3) = xs(3,1)*xs(1,2) - xs(1,1)*xs(3,2)
      xsinv(3,3) = xs(1,1)*xs(2,2) - xs(2,1)*xs(1,2)     

      do j = 1, 3
        do i = 1, 3
          xsinv(i,j) = xsinv(i,j)/xsj
        end do
      end do

!  Compute derivatives of shape functions with respect to nodal
!  coordinates.

!  shp(j,k) = partial N^k/partial r_j     (j=1,2,3)
!
!  Here it is calculated as

!              3   partial N^k   partial xi_i     3
!  shp(j,k) = Sum ------------- -------------- = Sum shpxi(i,k) xsinv(i,j)
!             i=1  partial xi_i  partial r_j     i=1

      do k = 1, 10
        do j = 1,3
          shp(j,k) = 0.0d0
          do i = 1,3
            shp(j,k) = shp(j,k)+shpxi(i,k)*xsinv(i,j)
          end do
        end do
      end do

!Check whether product of Jacobian and its inverse equals
!the unit matrix
!  do j = 1,3
!    do i = 1, 3
!       prod(i,j) = 0.0d0
!    end do
!  end do

!  do j = 1, 3
!    do i = 1, 3
!      do k = 1, 3
!         prod(i,j) = prod(i,j)+xs(k,i)*xsinv(k,j)
!      end do
!    end do
!  end do

!  write(iow,'('' product of Jacobian and its inverse ''/
! :     3(2X,E16.9,1X)/3(2X,E16.9,1X)/3(2X,E16.9,1X))')
! :     ((prod(i,j),j=1,3),i=1,3)           

      end
      
