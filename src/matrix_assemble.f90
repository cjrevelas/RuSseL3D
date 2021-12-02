subroutine matrix_assemble(Rg2_per_mon, wa)
!-----------------------------------------------------------------------------------------------------------!
use kcw
use geometry
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!   
integer                  :: lint
integer                  :: m_1, n_1
integer                  :: i, j, l, n, m, nn, i1
integer, dimension (nel) :: gl_index

real(8), intent(in)                   :: Rg2_per_mon
real(8), intent(in), dimension(numnp) :: wa
real(8)                               :: xsj
real(8), dimension(ndm,nel)           :: xl
real(8), dimension(4,11)              :: shp
real(8), dimension(5,11)              :: sv
!-----------------------------------------------------------------------------------------------------------!
i1=0

F_m%c=0.d0
F_m%k=0.d0
F_m%g=0.d0
F_m%w=0.d0

do nn = 1, numel
    !1. loop over all nodes of current element
    !2. find global index of node
    !3. copy coordinates from global array to local array concerning current element

    do i = 1, nel
        gl_index(i) = ix(i,nn)
        do j = 1, ndm
            xl(j,i) = xc(j, gl_index(i))
        enddo
    enddo

    !set up for gauss quadrature
    l=3
    call gausspoints(l, lint, sv)

    do l = 1, lint

        i1 = nel*nel*(nn-1)

        call tetshp(sv(1,l), xl, ndm, nel, xsj, shp)

        do m = 1,nel
            m_1 = gl_index(m)
            do n = 1,nel
                n_1 = gl_index(n)

                i1 = i1 + 1

                F_m%c(i1) = F_m%c(i1) + shp(4,n)*shp(4,m)*xsj*sv(5,l)

                F_m%k(i1) = F_m%k(i1) + Rg2_per_mon &
                                      * (shp(1,n)*shp(1,m)+shp(2,n)*shp(2,m)+shp(3,n)*shp(3,m))*xsj*sv(5,l)

                F_m%w(i1) = F_m%w(i1) + &
                                wa(n_1)*shp(4,n)*shp(4,m)*xsj*sv(5,l)
            enddo !n
        enddo !m
    enddo !l
enddo !nn

do i = 1, all_el
     if (con_l2(i)/=i) then
         F_m%k(con_l2(i)) = F_m%k(con_l2(i)) + F_m%k(i)
         F_m%k(i)=0.
         F_m%c(con_l2(i)) = F_m%c(con_l2(i)) + F_m%c(i)
         F_m%c(i)=0.
         F_m%w(con_l2(i)) = F_m%w(con_l2(i)) + F_m%w(i)
         F_m%w(i)=0.
     endif
enddo

#ifdef DEBUG_OUTPUTS
open(unit=400, file = 'matrix_assembly_kcw.out.txt')
write(400,'(3(A20))')"F_m%k","F_m%c","F_m%w"
do i = 1, all_el
    write(400,'(3(E20.9))')F_m%k(i), F_m%c(i), F_m%w(i)
enddo
close(400)
#endif

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine matrix_assemble
