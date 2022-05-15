!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine fem_bcs_and_nonzeros(ds, mumpsMatrixType, node_belongs_to_dirichlet_face)
!------------------------------------------------------------------------------------------------------!
use kcw_mod,         only: F_m, A_m, NNZ
use geometry_mod,    only: total_num_of_node_pairs, nel, numel, numnp, node_pairing_xx_hash, &
&                          node_pairing_yy_hash, node_pairing_zz_hash
use constants_mod,   only: tol
use parser_vars_mod, only: periodicAxisId
use flags_mod,       only: mumps_asymm, mumps_posDef, mumps_genSymm

!#define PRINT_AFULL
#ifdef PRINT_AFULL
use iofiles_mod, only: A_matrix_full
#endif
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: mumpsMatrixType
integer             :: ii, jj, kk

logical, intent(in), dimension(numnp) :: node_belongs_to_dirichlet_face
logical, dimension(nel*numel)         :: set_diag_to_one

real(8), intent(in) :: ds

#ifdef PRINT_AFULL
real(8), allocatable, dimension(:,:) :: A_full
#endif
!------------------------------------------------------------------------------------------------------!
F_m%g  = F_m%c + ds * (F_m%k + F_m%w)
F_m%rh = F_m%c

! Apply periodic boundary conditions
if (periodicAxisId(1)) call fem_apply_periodic_bcs(node_pairing_xx_hash)
if (periodicAxisId(2)) call fem_apply_periodic_bcs(node_pairing_yy_hash)
if (periodicAxisId(3)) call fem_apply_periodic_bcs(node_pairing_zz_hash)

! Apply Dirichlet boundary conditions
! TODO: needs generalization to nonzero values
! In case the matrix is symmetric, remove the zero lines and rows diagonal componets with Dirichlet BC q=0.
set_diag_to_one=.true.
if ((mumpsMatrixType.eq.mumps_posDef).or.(mumpsMatrixType.eq.mumps_genSymm)) then
    do kk = 1, total_num_of_node_pairs
        if (F_m%is_zero(kk)) cycle
        if (F_m%row(kk)==0)  cycle

        ii = F_m%row(kk)
        jj = F_m%col(kk)

        if (ii > jj) F_m%g(kk) = 0.0d0

        if (node_belongs_to_dirichlet_face(ii).or.node_belongs_to_dirichlet_face(jj)) then
            F_m%g(kk) = 0.0d0
            if (ii==jj.and.set_diag_to_one(ii)) then
                F_m%g(kk)           = 1.0d0
                set_diag_to_one(ii) = .false.
            endif
        endif
    enddo
endif

if (mumpsMatrixType.eq.mumps_asymm) then
    do kk = 1, total_num_of_node_pairs
        if (F_m%is_zero(kk)) cycle

        if (F_m%row(kk)==0) cycle

        jj = F_m%col(kk)
        ii = F_m%row(kk)

        if (node_belongs_to_dirichlet_face(ii)) then
            F_m%g(kk) = 0.0d0
            if (ii==jj.and.set_diag_to_one(ii)) then
                F_m%g(kk)           =  1.0d0
                set_diag_to_one(ii) = .false.
            endif
        endif
    enddo
endif

! Determine non_zero entries
NNZ = 0
do kk = 1, total_num_of_node_pairs
    if (ABS(F_m%g(kk)) > tol) NNZ = NNZ + 1
enddo

allocate(A_m%value(NNZ))
allocate(A_m%col(NNZ))
allocate(A_m%row(NNZ))

NNZ = 0
do kk = 1, total_num_of_node_pairs
    if (ABS(F_m%g(kk)) > tol) then
        NNZ = NNZ + 1

        A_m%value(NNZ) = F_m%g(kk)
        A_m%row(NNZ)   = F_m%row(kk)
        A_m%col(NNZ)   = F_m%col(kk)
    endif
enddo

#ifdef PRINT_AFULL
allocate(A_full(numnp,numnp))
A_full = 0.0d0

do kk = 1, NNZ
   ii = A_m%row(kk)
   jj = A_m%col(kk)

   A_full(ii,jj) = A_m%value(kk)
enddo

open(unit=255, file = A_matrix_full)
do ii = 1, numnp
   write(255,*) (A_full(ii,jj), jj = 1, numnp)
enddo
close(255)

deallocate(A_full)
#endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine fem_bcs_and_nonzeros
