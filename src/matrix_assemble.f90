subroutine matrix_assemble
!-----------------------------------------------------------------------------------------------------------!
use mdata
use xdata
use kcw
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!   
integer :: lint

integer :: m_1, n_1

integer :: i, j, l, n, m, nn !,ii

integer, dimension (nel) :: gl_index

real(8) :: xl(ndm,nel)!, ul(ndm,nel,7)

real(8) :: shp(4,11), sv(5,11), xsj
!-----------------------------------------------------------------------------------------------------------!
i1=0

do nn = 1, numel
    !1. Loop over all nodes of current element
    !2. Find global index of node
    !3. Copy coordinates from global array to local array concerning current element

    do i = 1, nel
        gl_index(i) = ix(i,nn)
        do j = 1, ndm
            xl(j,i) = xc(j, gl_index(i))
        enddo
    enddo

    !Set up for gauss quadrature
    l=3
    call gauss_3d(l, lint, sv)

    do l = 1, lint

        i1 = nel*nel*(nn-1)

        call tetshp( sv(1,l), xl, ndm, nel, xsj, shp )

        do m = 1,nel
            m_1 = gl_index(m)
            do n = 1,nel
                n_1 = gl_index(n)

                i1 = i1 + 1

                c_m%value(i1) = c_m%value(i1) + shp(4,n)*shp(4,m)*xsj*sv(5,l)

                k_m%value(i1) = k_m%value(i1) + &
                               (Rgyr**2.)*(shp(1,n)*shp(1,m)+shp(2,n)*shp(2,m)+shp(3,n)*shp(3,m))*xsj*sv(5,l)

                w_m%value(i1) = w_m%value(i1) + &
                                wa(n_1)*shp(4,n)*shp(4,m)*xsj*sv(5,l) 
            enddo !n
        enddo !m
    enddo !l
enddo !nn

do i = 1, all_el
     if (con_l2(i)/=i) then
         k_m%value(con_l2(i)) = k_m%value(con_l2(i)) + k_m%value(i)
         k_m%value(i)=0.
         c_m%value(con_l2(i)) = c_m%value(con_l2(i)) + c_m%value(i)
         c_m%value(i)=0.
         w_m%value(con_l2(i)) = w_m%value(con_l2(i)) + w_m%value(i)
         w_m%value(i)=0.
     endif
enddo

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine matrix_assemble

    
