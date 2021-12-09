subroutine dirichlet(ds, mumps_matrix_type, node_in_q0_face)
!------------------------------------------------------------------------------------------------------!
use kcw,       only: F_m, A_m, NNZ
use geometry,  only: all_el, nel, numel, numnp
use constants, only: tol
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: mumps_matrix_type
integer             :: jj, kk, ii

logical, intent(in), dimension(numnp) :: node_in_q0_face
logical, dimension(nel*numel)         :: set_diag_to_one

real(8), intent(in) :: ds

#ifdef PRINT_AFULL
real(8), allocatable, dimension(:,:) :: A_full
#endif
!------------------------------------------------------------------------------------------------------!
F_m%g  = F_m%c + ds * (F_m%k + F_m%w)
F_m%rh = F_m%c

set_diag_to_one=.true.

!in case the matrix is symmetric remove the zero lines and rows diagonal componets with Dirichlet BC q=0.
if (mumps_matrix_type.eq.1.or.mumps_matrix_type.eq.2) then
    do ii = 1, all_el
        kk = F_m%col(ii)
        jj = F_m%row(ii)
        if (jj > kk) then
            F_m%g(ii) = 0.d0
        endif
        if (node_in_q0_face(jj).or.node_in_q0_face(kk)) then
            F_m%g(ii) = 0.d0
            if (jj==kk.and.set_diag_to_one(jj)) then
                F_m%g(ii) = 1.d0
                set_diag_to_one(jj)=.false.
            endif
        endif
    enddo
endif

if (mumps_matrix_type.eq.0) then
    do ii = 1, all_el
        kk = F_m%col(ii)
        jj = F_m%row(ii)
        if (node_in_q0_face(jj)) then
            F_m%g(ii) = 0.d0
            if (jj==kk.and.set_diag_to_one(jj)) then
                F_m%g(ii) = 1.d0
                set_diag_to_one(jj)=.false.
            endif
        endif
    enddo
endif

!determine non_zero entries
NNZ = 0
do jj = 1, all_el
    if (ABS(F_m%g(jj)) > tol) then
        NNZ = NNZ + 1
    endif
enddo

allocate(A_m%value(NNZ))
allocate(A_m%col(NNZ))
allocate(A_m%row(NNZ))

NNZ = 0
do jj = 1, all_el
    if (ABS(F_m%g(jj)) > tol) then
        NNZ = NNZ + 1

        A_m%value(NNZ) = F_m%g(jj)
        A_m%row(NNZ)   = F_m%row(jj)
        A_m%col(NNZ)   = F_m%col(jj)
    endif
enddo

#ifdef PRINT_AFULL
allocate(A_full(numnp,numnp))
A_full = 0.d0

do ii = 1, NNZ
   kk = A_m%col(ii)
   jj = A_m%row(ii)
   A_full(jj,kk)=A_m%value(ii)
enddo

open(unit=255, file = A_matrix_full)
do jj = 1, numnp
   write(255,*)(A_full(jj,kk), kk = 1, numnp)
enddo
close(255)

deallocate(A_full)
#endif
!------------------------------------------------------------------------------------------------------!
end subroutine dirichlet
