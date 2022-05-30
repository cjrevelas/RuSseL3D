!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine fem_bcs_and_nonzeros(ds, mumpsMatrixType, nodeBelongsToDirichletFace)
!------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module

use kcw_mod,         only: F_m, A_m, NNZ
use geometry_mod,    only: numTotalNodePairs, numNodesLocalTypeDomain, numElementsTypeDomain, &
                           numNodes, nodePairingXXhash, nodePairingYYhash, nodePairingZZhash
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
integer             :: ii, jj, kk, mm, nn

logical, intent(in), dimension(numNodes)                          :: nodeBelongsToDirichletFace
logical, dimension(numNodesLocalTypeDomain*numElementsTypeDomain) :: set_diag_to_one
logical                                                           :: success

type(fhash_type_iterator__ints_double) :: nodePairingXXit, nodePairingXXit_aux, nodePairingYYit, nodePairingYYit_aux, nodePairingZZit
type(ints_type)                        :: nodePairingXXkey, nodePairingYYkey, nodePairingZZkey, dest_both_key
integer                                :: nodePairingXXvalue, nodePairingYYvalue, nodePairingZZvalue, dest_both

real(8), intent(in) :: ds

integer :: source_xx, source_yy, source_zz, source_aux=0
integer :: dest_xx, dest_yy, dest_zz, dest_triple=0

#ifdef PRINT_AFULL
real(8), allocatable, dimension(:,:) :: A_full
#endif
!------------------------------------------------------------------------------------------------------!
F_m%g  = F_m%c + ds * (F_m%k + F_m%w)
F_m%rh = F_m%c

! Apply periodic boundary conditions
if (periodicAxisId(1)) call fem_apply_periodic_bcs(nodePairingXXhash)
if (periodicAxisId(2)) call fem_apply_periodic_bcs(nodePairingYYhash)

if (periodicAxisId(1).AND.periodicAxisId(2)) then
  call nodePairingYYit%begin(nodePairingYYhash)

  allocate(dest_both_key%ints(1))

  do ii = 1, nodePairingYYhash%key_count()
    call nodePairingYYit%next(nodePairingYYkey, nodePairingYYvalue)
    source_yy = nodePairingYYkey%ints(1)  ! source_yy = 7
    dest_yy   = nodePairingYYvalue        ! dest_yy   = 1

    dest_both_key%ints(1) = dest_yy
    call nodePairingXXhash%get(dest_both_key, dest_both, success) ! dest_both = 38

    call nodePairingXXit%begin(nodePairingXXhash)
    do jj = 1, nodePairingXXhash%key_count()
      call nodePairingXXit%next(nodePairingXXkey, nodePairingXXvalue)
      source_xx = nodePairingXXkey%ints(1) ! source_xx = 7
      dest_xx   = nodePairingXXvalue       ! dest_xx   = 23

      if (source_yy.eq.source_xx) then
        do mm = 1, numTotalNodePairs
          if ((F_m%row(mm).eq.dest_xx).AND.(F_m%col(mm).eq.source_xx)) then ! (23,7)
            do nn = 1, numTotalNodePairs
              if ((F_m%row(nn).eq.dest_both).AND.(F_m%col(nn).eq.dest_yy)) then ! (38,1)
                F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
                F_m%g(nn) = 0.0d0
                exit
              endif
            enddo
          endif
        enddo
      endif
    enddo
  enddo
endif

if (periodicAxisId(3)) call fem_apply_periodic_bcs(nodePairingZZhash)

if (periodicAxisId(1).AND.periodicAxisId(2).AND.periodicAxisId(3)) then
  call nodePairingZZit%begin(nodePairingZZhash)

  do ii = 1, nodePairingZZhash%key_count()
    call nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)
    source_zz = nodePairingZZkey%ints(1) ! source_zz = 9 or 48
    dest_zz   = nodePairingZZvalue       ! dest_zz   = 1 or 38

    call nodePairingYYit%begin(nodePairingYYhash)
    do jj = 1, nodePairingYYhash%key_count()
      call nodePairingYYit%next(nodePairingYYkey, nodePairingYYvalue)
      source_yy = nodePairingYYkey%ints(1) ! source_yy = 7 or 23
      dest_yy   = nodePairingYYvalue       ! dest_yy   = 1 or 38

      if (dest_zz.eq.dest_yy) then
        call nodePairingYYit_aux%begin(nodePairingYYhash)
        do kk = 1, nodePairingYYhash%key_count()
          call  nodePairingYYit_aux%next(nodePairingYYkey, nodePairingYYvalue)
          if (nodePairingYYvalue.eq.source_zz) then
            source_aux = nodePairingYYkey%ints(1) ! source_aux = 21, 31
            exit
          endif
        enddo

        do mm = 1, numTotalNodePairs
          if ((F_m%row(mm).eq.source_zz).AND.(F_m%col(mm).eq.source_aux)) then
            do nn = 1, numTotalNodePairs
              if ((F_m%row(nn).eq.dest_zz).AND.(F_m%col(nn).eq.source_yy)) then
                F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
                F_m%g(nn) = 0.0d0
                exit
              endif
            enddo
          endif
        enddo

        call nodePairingXXit%begin(nodePairingXXhash)
        do kk = 1, nodePairingXXhash%key_count()
          call nodePairingXXit%next(nodePairingXXkey, nodePairingXXvalue)
          dest_xx = nodePairingXXvalue
          if (dest_zz.eq.dest_xx) dest_triple = dest_xx
        enddo
      endif
    enddo
  enddo

  call nodePairingZZit%begin(nodePairingZZhash)
  do ii = 1, nodePairingZZhash%key_count()
    call nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)
    source_zz = nodePairingZZkey%ints(1) ! source_zz = 31
    dest_zz   = nodePairingZZvalue       ! dest_zz   = 23

    call nodePairingXXit%begin(nodePairingXXhash)
    do jj = 1, nodePairingXXhash%key_count()
      call nodePairingXXit%next(nodePairingXXkey, nodePairingXXvalue)
      source_xx = nodePairingXXkey%ints(1) ! source_xx = 7
      dest_xx   = nodePairingXXvalue       ! dest_xx   = 23

      if ((dest_zz.eq.dest_xx).AND.(dest_zz.ne.dest_triple)) then
        !write(6,*) source_zz, dest_zz, source_xx, dest_xx
        call nodePairingXXit_aux%begin(nodePairingXXhash)
        do kk = 1, nodePairingXXhash%key_count()
          call nodePairingXXit_aux%next(nodePairingXXkey, nodePairingXXvalue)
          if (nodePairingXXvalue.eq.source_zz) then
            source_aux = nodePairingXXkey%ints(1) ! source_aux = 21
            exit
          endif
        enddo

        do mm = 1, numTotalNodePairs
          if ((F_m%row(mm).eq.source_zz).AND.(F_m%col(mm).eq.source_aux)) then
            do nn = 1, numTotalNodePairs
              if ((F_m%row(nn).eq.dest_zz).AND.(F_m%col(nn).eq.source_xx)) then
                F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
                F_m%g(nn) = 0.0d0
              endif
            enddo
          endif
        enddo
      endif
    enddo
  enddo
endif

! Prepare stiffness matrix for Dirichlet boundary conditions
! In case the matrix is symmetric, remove the zero lines and rows diagonal componets with Dirichlet BC q=0.
set_diag_to_one=.true.
if ((mumpsMatrixType.eq.mumps_posDef).or.(mumpsMatrixType.eq.mumps_genSymm)) then
  do kk = 1, numTotalNodePairs
    if (F_m%is_zero(kk)) cycle
    if (F_m%row(kk)==0)  cycle

    ii = F_m%row(kk)
    jj = F_m%col(kk)

    if (ii > jj) F_m%g(kk) = 0.0d0

    if (nodeBelongsToDirichletFace(ii).OR.nodeBelongsToDirichletFace(jj)) then
      F_m%g(kk) = 0.0d0
      if (ii==jj.and.set_diag_to_one(ii)) then
        F_m%g(kk)           = 1.0d0
        set_diag_to_one(ii) = .false.
      endif
    endif
  enddo
endif

if (mumpsMatrixType.eq.mumps_asymm) then
  do kk = 1, numTotalNodePairs
    if (F_m%is_zero(kk)) cycle

    if (F_m%row(kk)==0) cycle

    jj = F_m%col(kk)
    ii = F_m%row(kk)

    if (nodeBelongsToDirichletFace(ii)) then
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
do kk = 1, numTotalNodePairs
  if (ABS(F_m%g(kk)) > tol) NNZ = NNZ + 1
enddo

allocate(A_m%value(NNZ))
allocate(A_m%col(NNZ))
allocate(A_m%row(NNZ))

NNZ = 0
do kk = 1, numTotalNodePairs
  if (ABS(F_m%g(kk)) > tol) then
    NNZ = NNZ + 1

    A_m%value(NNZ) = F_m%g(kk)
    A_m%row(NNZ)   = F_m%row(kk)
    A_m%col(NNZ)   = F_m%col(kk)
  endif
enddo

#ifdef PRINT_AFULL
allocate(A_full(numNodes,numNodes))
A_full = 0.0d0

do kk = 1, NNZ
  ii = A_m%row(kk)
  jj = A_m%col(kk)

  A_full(ii,jj) = A_m%value(kk)
enddo

open(unit=255, file = A_matrix_full)
do ii = 1, numNodes
  write(255,*) (A_full(ii,jj), jj = 1, numNodes)
enddo
close(255)

deallocate(A_full)
#endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine fem_bcs_and_nonzeros
