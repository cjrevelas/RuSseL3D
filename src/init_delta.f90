!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_delta()
!-----------------------------------------------------------------------------------------------------------!
use geometry_mod,     only: numnp
use parser_vars_mod,  only: gr_exist
use arrays_mod,       only: phi_gr_indiv
use iofiles_mod,      only: graftFile
use delta_mod
use error_handing_mod
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: i1, iog
!-----------------------------------------------------------------------------------------------------------!
if (gr_exist.eq.1) then
    iog = 19

    open(unit=iog, file=graftFile)

    read(iog,*)
    read(iog,*)
    read(iog,*)
    read(iog,*) targetNumGraftedChains
    read(iog,*)
    read(iog,*)
    read(iog,*)
    read(iog,*)
    read(iog,*)

    allocate(gpid(targetNumGraftedChains), delta_numer(targetNumGraftedChains), gp_init_value(targetNumGraftedChains))

    gpid          = 0
    delta_numer   = 0.0d0
    gp_init_value = 0.0d0

    ! Specify grafting points
    do i1 = 1, targetNumGraftedChains
        read(iog,*) gpid(i1), gp_init_value(i1), delta_numer(i1)

        if (gpid(i1) > numnp) then
            write(ERROR_MESSAGE,'("ID of grafted chain (",I10,") is larger from numnp (",I10,"    )")') gpid(i1), numnp
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    enddo

    close(iog)
else
    targetNumGraftedChains = 0
    allocate(gpid(targetNumGraftedChains))
endif

allocate(phi_gr_indiv(numnp,targetNumGraftedChains))
!-----------------------------------------------------------------------------------------------------------!
end subroutine init_delta
