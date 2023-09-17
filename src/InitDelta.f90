!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine InitDelta()
!-----------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: numNodes
use parser_vars_mod,  only: graftedExist
use arrays_mod,       only: phiGraftedIndiv
use iofiles_mod,      only: IO_graftFile
use delta_mod
use error_handling_mod
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: graftedChain, iog
!-----------------------------------------------------------------------------------------------------------!
if (graftedExist.eq.1) then
  iog = 19

  open(unit=iog, file=IO_graftFile)
  read(iog,*)
  read(iog,*)
  read(iog,*)
  read(iog,*) targetNumGraftedChains
  read(iog,*)
  read(iog,*)
  read(iog,*)
  read(iog,*)
  read(iog,*)

  allocate(graftPointId(targetNumGraftedChains), deltaNumerical(targetNumGraftedChains), graftPointValue(targetNumGraftedChains))

  graftPointId    = 0
  deltaNumerical  = 0.0d0
  graftPointValue = 0.0d0

  do graftedChain = 1, targetNumGraftedChains
    read(iog,*) graftPointId(graftedChain), graftPointValue(graftedChain), deltaNumerical(graftedChain)

    if (graftPointId(graftedChain) > numNodes) then
      write(ERROR_MESSAGE,'("ID of grafted chain (",I10,") is larger from numNodes (",I10,"    )")') graftPointId(graftedChain), numNodes
      CALL ExitWithError(1, 1, 1, ERROR_MESSAGE)
    endif
  enddo

  close(iog)
else
  targetNumGraftedChains = 0
  allocate(graftPointId(targetNumGraftedChains))
endif

allocate(phiGraftedIndiv(numNodes,targetNumGraftedChains))
!-----------------------------------------------------------------------------------------------------------!
end subroutine InitDelta
