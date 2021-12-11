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
!----------------------------------------------------------------------------------------------------------------!
kk = 0

F_m%c = 0.0d0
F_m%k = 0.0d0
F_m%g = 0.0d0
F_m%w = 0.0d0

! Assembly element matrices
do elem = 1, numel
    ! 1. Loop over all nodes of current element
    ! 2. Find global index of node
    ! 3. Copy coordinates from global array to local array concerning current element

    do ii = 1, nel
        global_index(ii) = global_node_id_type_domain(ii,elem)
        do jj = 1, ndm
            xl(jj,ii) = xc(jj, global_index(ii))
        enddo
    enddo

    ! Set up for gauss quadrature
    ll = 3
    call fem_gausspoints(ll, lint, sv)

    do ll = 1, lint

        kk = nel * nel * (elem-1)     ! This index goes from zero to num_of_bulk_pairs = nel * nel * numel

        call fem_tetshpfun(sv(1,ll), xl, ndm, nel, xsj, shp)

        do mm = 1, nel
            do nn = 1, nel
                pp = global_index(nn)

                kk = kk + 1

                F_m%c(kk) = F_m%c(kk) + shp(4,nn)*shp(4,mm)*xsj*sv(5,ll)

                F_m%k(kk) = F_m%k(kk) + Rg2_per_mon &
                                      * (shp(1,nn)*shp(1,mm)+shp(2,nn)*shp(2,mm)+shp(3,nn)*shp(3,mm))*xsj*sv(5,ll)

                F_m%w(kk) = F_m%w(kk) + wa(pp)*shp(4,nn)*shp(4,mm)*xsj*sv(5,ll)
            enddo !nn
        enddo !mm
    enddo !ll
enddo !elem

! Assembly global matrix using element matrices and node_pair_id hash matrix created in parser_mesh.f90
! TODO: Lots of nonzero entries are created with the following process -> consider removal
do kk = 1, num_of_bulk_pairs
    if (F_m%is_zero(kk)) then
        ! Add up contributions of same pairs met multiple times
        F_m%k(node_pair_id(kk)) = F_m%k(node_pair_id(kk)) + F_m%k(kk)
        F_m%k(kk)               = 0.0d0
        F_m%c(node_pair_id(kk)) = F_m%c(node_pair_id(kk)) + F_m%c(kk)
        F_m%c(kk)               = 0.0d0
        F_m%w(node_pair_id(kk)) = F_m%w(node_pair_id(kk)) + F_m%w(kk)
        F_m%w(kk)               = 0.0d0
    endif
enddo

#ifdef DEBUG_OUTPUTS
open(unit=400, file = matrix_assembly)
write(400,'(3(2X,A16))') "F_m%k","F_m%c","F_m%w"
do kk = 1, num_of_bulk_pairs
    write(400,'(3(2X,E16.9))')F_m%k(kk), F_m%c(kk), F_m%w(kk)
enddo
close(400)
#endif

return
!----------------------------------------------------------------------------------------------------------------!
end subroutine fem_matrix_assemble
