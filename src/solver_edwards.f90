!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine solver_edwards(ds, ns, mumpsMatrixType, q, q_final, nodeBelongsToDirichletFace)
!----------------------------------------------------------------------------------------------------------!
use kcw_mod,         only: A_m, F_m, rdiag1
use parser_vars_mod, only: numDirichletFaces, numNanoparticleFaces,   &
                           dirichletFaceValue, nanoparticleFaceValue, &
                           dirichletFaceId, dirichletFaceValue,       &
                           nanoparticleFaceId, nanoparticleFaceValue
use geometry_mod,    only: numNodes, numTotalNodePairs, nodeBelongsToFaceId
#ifdef USE_MPI
use mpistuff
#endif
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
#ifdef USE_MPI
include "mpif.h"
#endif
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns, mumpsMatrixType
integer             :: ii, jj, kk, time_step, tools_sys_time, t_init, t_final, face

logical, intent(in), dimension(numNodes) :: nodeBelongsToDirichletFace

real(8), intent(in), dimension(ns+1)          :: ds
real(8), intent(inout), dimension(2,numNodes)    :: q
real(8), intent(inout), dimension(ns+1,numNodes) :: q_final
!----------------------------------------------------------------------------------------------------------!
t_init = tools_sys_time()

do time_step = 2, ns+1
  call fem_bcs_and_nonzeros(ds(time_step), mumpsMatrixType, nodeBelongsToDirichletFace)

#ifdef USE_MPI
  ! Send a continue (.true.) signal to the slaves
  call MPI_BCAST(.True., 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

  if (ns.ge.10.and.MOD(time_step+1,ns/10).eq.0) then
    write(6,'(I3,1X)',advance='no') NINT((time_step-2.0d0)/ns*100.0d0)
  elseif (ns.lt.10.and.MOD(time_step+1,1).eq.0) then
    write(6,'(I3,1X)',advance='no') NINT((time_step-2.0d0)/ns*100.0d0)
  endif

  ! Form the RHS of the linear system of equations to be solved
  rdiag1 = 0.0d0

  do kk = 1, numTotalNodePairs
    if (F_m%is_zero(kk)) cycle
    if (F_m%row(kk)==0) cycle ! TEMP SOLUTION

    ii = F_m%row(kk)
    jj = F_m%col(kk)

    rdiag1(ii) = rdiag1(ii) + F_m%rh(kk)*q(1,jj)
  enddo

  ! Assing value at the Dirichlet boundaries
  do face = 1, numDirichletFaces
    do ii = 1, numNodes
      if (nodeBelongsToFaceId(ii) == dirichletFaceId(face)) rdiag1(ii) = dirichletFaceValue(face)
    enddo
  enddo

  do face = 1, numNanoparticleFaces
    do ii = 1, numNodes
      if (nodeBelongsToFaceId(ii) == nanoparticleFaceId(face)) rdiag1(ii) = nanoparticleFaceValue(face)
    enddo
  enddo

  ! Solve the linear system of equations
  call solver_mumps(mumpsMatrixType)

  ! Update solution/propagator
  do kk = 1,numNodes
    q(2,kk) = rdiag1(kk)
  enddo

  ! Save propagators for convolution
  do kk = 1,numNodes
    q_final(time_step,kk) = q(2,kk)
    q(1,kk)               = q(2,kk)
  enddo

  ! Deallocate the A arrays generated by dirichlet subroutine
  deallocate(A_m%value)
  deallocate(A_m%col)
  deallocate(A_m%row)
enddo

write(6,'(1X,I3)',advance='no') 100

t_final = tools_sys_time()

write(6,'(" :",I6)') t_final - t_init

return
!----------------------------------------------------------------------------------------------------------!
end subroutine solver_edwards
