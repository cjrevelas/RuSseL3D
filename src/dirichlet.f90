subroutine dirichlet(ds, ns, mumps_matrix_type)
!------------------------------------------------------------------------------------------------------!
use kcw
use geometry
use constants
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns, mumps_matrix_type
integer             :: i, f, i1

logical, dimension(nel*numel) :: set_diag_to_one

real(8), intent(in), dimension(ns+1) :: ds

#ifdef PRINT_AFULL
real(8), allocatable, dimension(:,:) :: A_full
#endif
!------------------------------------------------------------------------------------------------------!
F_m%g  = F_m%c + ds(1)*(F_m%k + F_m%w)
F_m%rh = F_m%c

set_diag_to_one=.true.

!in case the matrix is symmetric remove the zero lines and rows
!diagonal componets with Dirichlet BC q=0.
if (mumps_matrix_type.eq.1.or.mumps_matrix_type.eq.2) then
    do i1 = 1, all_el
        f = F_m%col(i1)
        i = F_m%row(i1)
        if (i > f) then
            F_m%g(i1) = 0.d0
        endif
        if (elem_in_q0_face(i).or.elem_in_q0_face(f)) then
            F_m%g(i1) = 0.d0
            if (i==f.and.set_diag_to_one(i)) then
                F_m%g(i1) = 1.d0
                set_diag_to_one(i)=.false.
            endif
        endif
    enddo
endif

if (mumps_matrix_type.eq.0) then
    do i1 = 1, all_el
        f = F_m%col(i1)
        i = F_m%row(i1)
        if (elem_in_q0_face(i)) then
            F_m%g(i1) = 0.d0
            if (i==f.and.set_diag_to_one(i)) then
                F_m%g(i1) = 1.d0
                set_diag_to_one(i)=.false.
            endif
        endif
    enddo
endif
!***********************DETERMINE NON-ZERO ENTRIES*******************!
NNZ = 0
do i = 1, all_el
    if (abs(F_m%g(i)) > tol) then
        NNZ = NNZ + 1
    endif
enddo

allocate(A_m%value(NNZ))
allocate(A_m%col(NNZ))
allocate(A_m%row(NNZ))

NNZ = 0
do i = 1, all_el
    if (abs(F_m%g(i)) > tol) then
        NNZ = NNZ + 1

        A_m%value(NNZ) = F_m%g(i)
        A_m%row(NNZ)   = F_m%row(i)
        A_m%col(NNZ)   = F_m%col(i)
    endif
enddo

#ifdef PRINT_AFULL
allocate(A_full(numnp,numnp))
A_full = 0.d0
do i1 = 1, NNZ
   f = A_m%col(i1)
   i = A_m%row(i1)
   A_full(i,f)=A_m%value(i1)
enddo
open(unit=255, file = 'A_full.out.txt')
do i = 1, numnp
   write(255,*)(A_full(i,f), f = 1, numnp)
enddo
close(255)
deallocate(A_full)
#endif
!------------------------------------------------------------------------------------------------------!
end subroutine dirichlet
