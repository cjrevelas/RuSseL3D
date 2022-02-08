!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine fem_bcs_and_nonzeros(ds, mumps_matrix_type, node_belongs_to_dirichlet_face)
!------------------------------------------------------------------------------------------------------!
use kcw_mod,         only: F_m, A_m, NNZ
use geometry_mod,    only: total_num_of_node_pairs, nel, numel, numnp,                          &
&                      node_pairing_xx_hash, node_pairing_xx_it, node_pairing_xx_key, node_pairing_xx_value, &
&                      node_pairing_yy_hash, node_pairing_yy_it, node_pairing_yy_key, node_pairing_yy_value, &
&                      node_pairing_zz_hash, node_pairing_zz_it, node_pairing_zz_key, node_pairing_zz_value
use constants_mod,   only: tol
use parser_vars_mod, only: periodic_axis_id

!#define PRINT_AFULL
#ifdef PRINT_AFULL
use iofiles_mod, only: A_matrix_full
#endif
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: mumps_matrix_type
integer             :: ii, jj, kk, mm, nn
integer             :: source_xx, dest_xx
integer             :: source_yy, dest_yy
integer             :: source_zz, dest_zz

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
! xx
if (periodic_axis_id(1)) then
    call node_pairing_xx_it%begin(node_pairing_xx_hash)

    do kk = 1, node_pairing_xx_hash%key_count()
        call node_pairing_xx_it%next(node_pairing_xx_key, node_pairing_xx_value)

        source_xx = node_pairing_xx_key%ints(1)
        dest_xx   = node_pairing_xx_value

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%col(mm).eq.source_xx).and.(F_m%row(mm).ne.dest_xx)) then
                do nn = 1, total_num_of_node_pairs
                    if (F_m%is_zero(nn)) cycle
                    if (F_m%row(nn)==0)  cycle

                    if ((F_m%row(nn).eq.F_m%row(mm)).and.(F_m%col(nn).eq.dest_xx)) then
                        F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                        F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                    endif
                enddo
            endif

            if ((F_m%row(mm).eq.source_xx).and.(F_m%col(mm).ne.dest_xx)) then
                do nn = 1, total_num_of_node_pairs
                    if (F_m%is_zero(nn)) cycle
                    if (F_m%row(nn)==0)  cycle

                    if ((F_m%row(nn).eq.dest_xx).and.(F_m%col(nn).eq.F_m%col(mm))) then
                        F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                        F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                    endif
                enddo
            endif
        enddo

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%row(mm).eq.source_xx).and.(F_m%col(mm).eq.source_xx)) then
                do nn = 1, total_num_of_node_pairs
                    if (F_m%is_zero(nn)) cycle
                    if (F_m%row(nn)==0)  cycle

                    if ((F_m%row(nn).eq.dest_xx).and.(F_m%col(nn).eq.dest_xx)) then
                        F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                        F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                    endif
                enddo
            endif
        enddo

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%row(mm).eq.dest_xx).or.(F_m%col(mm).eq.dest_xx)) then
                F_m%g(mm)  = 0.0
                F_m%rh(mm) = 0.0
            endif
        enddo

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%row(mm).eq.dest_xx).and.(F_m%col(mm).eq.dest_xx))   F_m%g(mm) =  1.0
            if ((F_m%row(mm).eq.dest_xx).and.(F_m%col(mm).eq.source_xx)) F_m%g(mm) = -1.d0
        enddo
    enddo
endif

! yy
if (periodic_axis_id(2)) then
    call node_pairing_yy_it%begin(node_pairing_yy_hash)

    do kk = 1, node_pairing_yy_hash%key_count()
        call node_pairing_yy_it%next(node_pairing_yy_key, node_pairing_yy_value)

        source_yy = node_pairing_yy_key%ints(1)
        dest_yy   = node_pairing_yy_value

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%col(mm).eq.source_yy).and.(F_m%row(mm).ne.dest_yy)) then
                do nn = 1, total_num_of_node_pairs
                    if (F_m%is_zero(nn)) cycle
                    if (F_m%row(nn)==0)  cycle

                    if ((F_m%row(nn).eq.F_m%row(mm)).and.(F_m%col(nn).eq.dest_yy)) then
                        F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                        F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                    endif
                enddo
            endif

            if ((F_m%row(mm).eq.source_yy).and.(F_m%col(mm).ne.dest_yy)) then
                do nn = 1, total_num_of_node_pairs
                    if (F_m%is_zero(nn)) cycle
                    if (F_m%row(nn)==0)  cycle

                    if ((F_m%row(nn).eq.dest_yy).and.(F_m%col(nn).eq.F_m%col(mm))) then
                        F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                        F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                    endif
                enddo
            endif
        enddo

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%row(mm).eq.source_yy).and.(F_m%col(mm).eq.source_yy)) then
                do nn = 1, total_num_of_node_pairs
                    if (F_m%is_zero(nn)) cycle
                    if (F_m%row(nn)==0)  cycle

                    if ((F_m%row(nn).eq.dest_yy).and.(F_m%col(nn).eq.dest_yy)) then
                        F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                        F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                    endif
                enddo
            endif
        enddo

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%row(mm).eq.dest_yy).or.(F_m%col(mm).eq.dest_yy)) then
                F_m%g(mm)  = 0.0
                F_m%rh(mm) = 0.0
            endif
        enddo

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%row(mm).eq.dest_yy).and.(F_m%col(mm).eq.dest_yy))   F_m%g(mm) =  1.0
            if ((F_m%row(mm).eq.dest_yy).and.(F_m%col(mm).eq.source_yy)) F_m%g(mm) = -1.d0
        enddo
    enddo
endif

! zz
if (periodic_axis_id(3)) then
    call node_pairing_zz_it%begin(node_pairing_zz_hash)

    do kk = 1, node_pairing_zz_hash%key_count()
        call node_pairing_zz_it%next(node_pairing_zz_key, node_pairing_zz_value)

        source_zz = node_pairing_zz_key%ints(1)
        dest_zz   = node_pairing_zz_value

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%col(mm).eq.source_zz).and.(F_m%row(mm).ne.dest_zz)) then
                do nn = 1, total_num_of_node_pairs
                    if (F_m%is_zero(nn)) cycle
                    if (F_m%row(nn)==0)  cycle

                    if ((F_m%row(nn).eq.F_m%row(mm)).and.(F_m%col(nn).eq.dest_zz)) then
                        F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                        F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                    endif
                enddo
            endif

            if ((F_m%row(mm).eq.source_zz).and.(F_m%col(mm).ne.dest_zz)) then
                do nn = 1, total_num_of_node_pairs
                    if (F_m%is_zero(nn)) cycle
                    if (F_m%row(nn)==0)  cycle

                    if ((F_m%row(nn).eq.dest_zz).and.(F_m%col(nn).eq.F_m%col(mm))) then
                        F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                        F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                    endif
                enddo
            endif
        enddo

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%row(mm).eq.source_zz).and.(F_m%col(mm).eq.source_zz)) then
                do nn = 1, total_num_of_node_pairs
                    if (F_m%is_zero(nn)) cycle
                    if (F_m%row(nn)==0)  cycle

                    if ((F_m%row(nn).eq.dest_zz).and.(F_m%col(nn).eq.dest_zz)) then
                        F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                        F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                    endif
                enddo
            endif
        enddo

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%row(mm).eq.dest_zz).or.(F_m%col(mm).eq.dest_zz)) then
                F_m%g(mm)  = 0.0
                F_m%rh(mm) = 0.0
            endif
        enddo

        do mm = 1, total_num_of_node_pairs
            if (F_m%is_zero(mm)) cycle
            if (F_m%row(mm)==0)  cycle

            if ((F_m%row(mm).eq.dest_zz).and.(F_m%col(mm).eq.dest_zz))   F_m%g(mm) =  1.0
            if ((F_m%row(mm).eq.dest_zz).and.(F_m%col(mm).eq.source_zz)) F_m%g(mm) = -1.d0
        enddo
    enddo
endif

set_diag_to_one=.true.
! In case the matrix is symmetric, remove the zero lines and rows diagonal componets with Dirichlet BC q=0.
if ((mumps_matrix_type.eq.1).or.(mumps_matrix_type.eq.2)) then
    do kk = 1, total_num_of_node_pairs
        if (F_m%is_zero(kk)) cycle
        if (F_m%row(kk)==0)  cycle

        ii = F_m%row(kk)
        jj = F_m%col(kk)

        if (ii > jj) F_m%g(kk) = 0.d0

        if (node_belongs_to_dirichlet_face(ii).or.node_belongs_to_dirichlet_face(jj)) then
            F_m%g(kk) = 0.d0
            if (ii==jj.and.set_diag_to_one(ii)) then
                F_m%g(kk)           = 1.d0
                set_diag_to_one(ii) = .false.
            endif
        endif
    enddo
endif

if (mumps_matrix_type.eq.0) then
    do kk = 1, total_num_of_node_pairs
        if (F_m%is_zero(kk)) cycle

        if (F_m%row(kk)==0) cycle

        jj = F_m%col(kk)
        ii = F_m%row(kk)

        if (node_belongs_to_dirichlet_face(ii)) then
            F_m%g(kk) = 0.d0
            if (ii==jj.and.set_diag_to_one(ii)) then
                F_m%g(kk)           =  1.d0
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
A_full = 0.d0

do kk = 1, NNZ
   ii = A_m%row(kk)
   jj = A_m%col(kk)

   A_full(ii,jj) = A_m%value(kk)
enddo

open(unit=255, file = A_matrix_full)
do ii = 1, numnp
   write(255,*)(A_full(ii,jj), jj = 1, numnp)
enddo
close(255)

deallocate(A_full)
#endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine fem_bcs_and_nonzeros
