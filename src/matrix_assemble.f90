subroutine matrix_assemble(Rg2_per_mon, wa)
!-----------------------------------------------------------------------------------------------------------!
use kcw,      only: F_m
use geometry, only: numnp, numel, ndm, nel, all_el, con_l2, ix, xc
use iofiles,  only: matrix_assembly
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer                  :: lint
integer                  :: i, j, l, n, m, nn, ii, kk
integer, dimension (nel) :: gl_index

real(8), intent(in)                   :: Rg2_per_mon
real(8), intent(in), dimension(numnp) :: wa
real(8)                               :: xsj
real(8), dimension(ndm,nel)           :: xl
real(8), dimension(4,11)              :: shp
real(8), dimension(5,11)              :: sv
!-----------------------------------------------------------------------------------------------------------!
ii = 0

F_m%c = 0.d0
F_m%k = 0.d0
F_m%g = 0.d0
F_m%w = 0.d0

!assembly element matrices
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

        ii = nel*nel*(nn-1)     !this index goes from zero to all_el=nel*nel*numel

        call tetshp(sv(1,l), xl, ndm, nel, xsj, shp)

        do m = 1,nel
        do n = 1,nel
            kk = gl_index(n)

            ii = ii + 1

            F_m%c(ii) = F_m%c(ii) + shp(4,n)*shp(4,m)*xsj*sv(5,l)

            F_m%k(ii) = F_m%k(ii) + Rg2_per_mon &
                                  * (shp(1,n)*shp(1,m)+shp(2,n)*shp(2,m)+shp(3,n)*shp(3,m))*xsj*sv(5,l)

            F_m%w(ii) = F_m%w(ii) + wa(kk)*shp(4,n)*shp(4,m)*xsj*sv(5,l)
        enddo !n
        enddo !m
    enddo !l
enddo !nn

!assembly global matrix using element matrices and con_12 hash matrix created in mesh.f90
do i = 1, all_el
    if (con_l2(i)/=i) then
        !add up contributions of same pairs met multiple times
        F_m%k(con_l2(i)) = F_m%k(con_l2(i)) + F_m%k(i)
        F_m%k(i)=0.
        F_m%c(con_l2(i)) = F_m%c(con_l2(i)) + F_m%c(i)
        F_m%c(i)=0.
        F_m%w(con_l2(i)) = F_m%w(con_l2(i)) + F_m%w(i)
        F_m%w(i)=0.
    endif
enddo

#ifdef DEBUG_OUTPUTS
open(unit=400, file = matrix_assembly)
write(400,'(A15,2A20)') "F_m%k","F_m%c","F_m%w"
do i = 1, all_el
    write(400,'(3(E20.9))')F_m%k(i), F_m%c(i), F_m%w(i)
enddo
close(400)
#endif

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine matrix_assemble
