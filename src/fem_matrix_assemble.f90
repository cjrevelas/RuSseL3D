!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine fem_matrix_assemble(Rg2_per_mon, wa)
!----------------------------------------------------------------------------------------------------------------!
use kcw_mod,      only: F_m
use geometry_mod, only: numnp, numel, ndm, nel, num_of_bulk_pairs, node_pair_id, global_node_id_type_domain, xc
use iofiles_mod,  only: matrix_assembly
!----------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------!
integer                 :: lint, elem
integer                 :: ii, jj, kk, ll, mm, nn, pp
integer, dimension(nel) :: global_index

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

! Assembly element matrices
do nn = 1, numel
    ! 1. Loop over all nodes of current element
    ! 2. Find global index of node
    ! 3. Copy coordinates from global array to local array concerning current element

    do i = 1, nel
        global_index(i) = global_node_id_type_domain(i,nn)
        do j = 1, ndm
            xl(j,i) = xc(j, global_index(i))
        enddo
    enddo

    ! Set up for gauss quadrature
    l=3
    call fem_gausspoints(l, lint, sv)

    do l = 1, lint

        ii = nel*nel*(nn-1)     ! This index goes from zero to all_el=nel*nel*numel

        call fem_tetshpfun(sv(1,l), xl, ndm, nel, xsj, shp)

        do m = 1, nel
            do n = 1, nel
                kk = global_index(n)

                ii = ii + 1

                F_m%c(ii) = F_m%c(ii) + shp(4,n)*shp(4,m)*xsj*sv(5,l)

                F_m%k(ii) = F_m%k(ii) + Rg2_per_mon &
                                      * (shp(1,n)*shp(1,m)+shp(2,n)*shp(2,m)+shp(3,n)*shp(3,m))*xsj*sv(5,l)

                F_m%w(ii) = F_m%w(ii) + wa(kk)*shp(4,n)*shp(4,m)*xsj*sv(5,l)
            enddo !n
        enddo !m
    enddo !l
enddo !nn

! Assembly global matrix using element matrices and con_12 hash matrix created in parser_mesh.f90
do i = 1, all_el
    if (con_l2(i)/=i) then
        ! Add up contributions of same pairs met multiple times
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
write(400,'(3(2X,A16))') "F_m%k","F_m%c","F_m%w"
do i = 1, all_el
    write(400,'(3(2X,E16.9))')F_m%k(i), F_m%c(i), F_m%w(i)
enddo
close(400)
#endif

return
!-----------------------------------------------------------------------------------------------------------!
end subroutine fem_matrix_assemble
