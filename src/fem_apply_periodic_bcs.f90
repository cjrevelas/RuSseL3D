!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine fem_apply_periodic_bcs(node_pairing_hash)
!----------------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module

use kcw_mod, only: F_m
use geometry_mod, only: total_num_of_node_pairs
!----------------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------------!
type(fhash_type__ints_double), intent(inout) :: node_pairing_hash
type(fhash_type_iterator__ints_double)       :: node_pairing_it
type(ints_type)                              :: node_pairing_key
integer                                      :: node_pairing_value

integer :: kk, mm, nn, source, dest
!----------------------------------------------------------------------------------------------------------------------------------!
call node_pairing_it%begin(node_pairing_hash)

! TODO: optimize periodic boundary conditions

do kk = 1, node_pairing_hash%key_count()
    call node_pairing_it%next(node_pairing_key, node_pairing_value)

    source = node_pairing_key%ints(1)
    dest   = node_pairing_value

    do mm = 1, total_num_of_node_pairs
        if (F_m%is_zero(mm)) cycle
        if (F_m%row(mm)==0)  cycle

        if ((F_m%col(mm).eq.source).and.(F_m%row(mm).ne.dest)) then
            do nn = 1, total_num_of_node_pairs
                if (F_m%is_zero(nn)) cycle
                if (F_m%row(nn)==0)  cycle

                if ((F_m%row(nn).eq.F_m%row(mm)).and.(F_m%col(nn).eq.dest)) then
                    F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                    F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                endif
            enddo
        endif

        if ((F_m%row(mm).eq.source).and.(F_m%col(mm).ne.dest)) then
            do nn = 1, total_num_of_node_pairs
                if (F_m%is_zero(nn)) cycle
                if (F_m%row(nn)==0)  cycle

                if ((F_m%row(nn).eq.dest).and.(F_m%col(nn).eq.F_m%col(mm))) then
                    F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                    F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                endif
            enddo
        endif
    enddo

    do mm = 1, total_num_of_node_pairs
        if (F_m%is_zero(mm)) cycle
        if (F_m%row(mm)==0)  cycle

        if ((F_m%row(mm).eq.source).and.(F_m%col(mm).eq.source)) then
            do nn = 1, total_num_of_node_pairs
                if (F_m%is_zero(nn)) cycle
                if (F_m%row(nn)==0)  cycle

                if ((F_m%row(nn).eq.dest).and.(F_m%col(nn).eq.dest)) then
                    F_m%g(mm)  = F_m%g(mm)  + F_m%g(nn)
                    F_m%rh(mm) = F_m%rh(mm) + F_m%rh(nn)
                endif
            enddo
        endif
    enddo

    do mm = 1, total_num_of_node_pairs
        if (F_m%is_zero(mm)) cycle
        if (F_m%row(mm)==0)  cycle

        if ((F_m%row(mm).eq.dest).or.(F_m%col(mm).eq.dest)) then
            F_m%g(mm)  = 0.0d0
            F_m%rh(mm) = 0.0d0
        endif

        if ((F_m%row(mm).eq.dest).and.(F_m%col(mm).eq.dest))   F_m%g(mm) =  1.0d0
        if ((F_m%row(mm).eq.dest).and.(F_m%col(mm).eq.source)) F_m%g(mm) = -1.0d0
    enddo
enddo

return
!----------------------------------------------------------------------------------------------------------------------------------!
end subroutine fem_apply_periodic_bcs
