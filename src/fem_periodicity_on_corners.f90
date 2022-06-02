subroutine fem_periodicity_on_corners()
!------------------------------------------------------------------------------------------------------!
use iso_fortran_env
use fhash_module__ints_double
use ints_module
use kcw_mod,      only: F_m
use geometry_mod, only: nodePairingXXhash, &
                        nodePairingYYhash, &
                        nodePairingZZhash, numTotalNodePairs
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
type(fhash_type_iterator__ints_double) :: nodePairingXXit, nodePairingYYit, nodePairingZZit
type(fhash_type_iterator__ints_double) :: nodePairingXXitAux, nodePairingYYitAux
type(ints_type)                        :: nodePairingXXkey, nodePairingYYkey, nodePairingZZkey
integer                                :: nodePairingXXvalue, nodePairingYYvalue, nodePairingZZvalue

integer :: sourceXX, destXX
integer :: sourceYY, destYY
integer :: sourceZZ, destZZ
integer :: sourceAux=0, destTriple=0
integer :: ii, jj, kk, mm, nn
!------------------------------------------------------------------------------------------------------!
call nodePairingZZit%begin(nodePairingZZhash)

do ii = 1, nodePairingZZhash%key_count()
  call nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)
  sourceZZ = nodePairingZZkey%ints(1) ! sourceZZ = 9 or 48
  destZZ   = nodePairingZZvalue       ! destZZ   = 1 or 38

  call nodePairingYYit%begin(nodePairingYYhash)
  do jj = 1, nodePairingYYhash%key_count()
    call nodePairingYYit%next(nodePairingYYkey, nodePairingYYvalue)
    sourceYY = nodePairingYYkey%ints(1) ! sourceYY = 7 or 23
    destYY   = nodePairingYYvalue       ! destYY   = 1 or 38

    if (destZZ.eq.destYY) then
      call nodePairingYYitAux%begin(nodePairingYYhash)
      do kk = 1, nodePairingYYhash%key_count()
        call  nodePairingYYitAux%next(nodePairingYYkey, nodePairingYYvalue)
        if (nodePairingYYvalue.eq.sourceZZ) then
          sourceAux = nodePairingYYkey%ints(1) ! sourceAux = 21, 31
          exit
        endif
      enddo

      do mm = 1, numTotalNodePairs
        if ((F_m%row(mm).eq.sourceZZ).AND.(F_m%col(mm).eq.sourceAux)) then
          do nn = 1, numTotalNodePairs
            if ((F_m%row(nn).eq.destZZ).AND.(F_m%col(nn).eq.sourceYY)) then
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
        destXX = nodePairingXXvalue
        if (destZZ.eq.destXX) destTriple = destXX
      enddo
    endif
  enddo
enddo

call nodePairingZZit%begin(nodePairingZZhash)
do ii = 1, nodePairingZZhash%key_count()
  call nodePairingZZit%next(nodePairingZZkey, nodePairingZZvalue)
  sourceZZ = nodePairingZZkey%ints(1) ! sourceZZ = 31
  destZZ   = nodePairingZZvalue       ! destZZ   = 23

  call nodePairingXXit%begin(nodePairingXXhash)
  do jj = 1, nodePairingXXhash%key_count()
    call nodePairingXXit%next(nodePairingXXkey, nodePairingXXvalue)
    sourceXX = nodePairingXXkey%ints(1) ! sourceXX = 7
    destXX   = nodePairingXXvalue       ! destXX   = 23

    if ((destZZ.eq.destXX).AND.(destZZ.ne.destTriple)) then
      call nodePairingXXitAux%begin(nodePairingXXhash)
      do kk = 1, nodePairingXXhash%key_count()
        call nodePairingXXitAux%next(nodePairingXXkey, nodePairingXXvalue)
        if (nodePairingXXvalue.eq.sourceZZ) then
          sourceAux = nodePairingXXkey%ints(1) ! sourceAux = 21
          exit
        endif
      enddo

      do mm = 1, numTotalNodePairs
        if ((F_m%row(mm).eq.sourceZZ).AND.(F_m%col(mm).eq.sourceAux)) then
          do nn = 1, numTotalNodePairs
            if ((F_m%row(nn).eq.destZZ).AND.(F_m%col(nn).eq.sourceXX)) then
              F_m%g(mm) = F_m%g(mm) + F_m%g(nn)
              F_m%g(nn) = 0.0d0
            endif
          enddo
        endif
      enddo
    endif
  enddo
enddo

return
!------------------------------------------------------------------------------------------------------!
end subroutine fem_periodicity_on_corners
