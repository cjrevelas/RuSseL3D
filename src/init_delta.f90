!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_delta()
!-----------------------------------------------------------------------------------------------------------!
use geometry,     only: numnp
use parser_vars,  only: gr_exist
use arrays,       only: phia_gr_indiv
use iofiles,      only: gp_filename
use delta
use error_handing
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: i1, iog
!-----------------------------------------------------------------------------------------------------------!
if (gr_exist.eq.1) then
    iog = 19

    open(unit=iog, file=gp_filename)

    read(iog,*)
    read(iog,*)
    read(iog,*)
    read(iog,*) num_gpoints
    read(iog,*)
    read(iog,*)
    read(iog,*)
    read(iog,*)
    read(iog,*)

    allocate(gpid(num_gpoints), delta_numer(num_gpoints), gp_init_value(num_gpoints))

    gpid          = 0
    delta_numer   = 0.d0
    gp_init_value = 0.d0

    ! Specify grafting points
    do i1 = 1, num_gpoints
        read(iog,*) gpid(i1), gp_init_value(i1), delta_numer(i1)

        if (gpid(i1) > numnp) then
            write(ERROR_MESSAGE,'("ID of grafted chain (",I10,") is larger from numnp (",I10,"    )")') gpid(i1), numnp
            call exit_with_error(1,1,1,ERROR_MESSAGE)
        endif
    enddo

    close(iog)
else
    num_gpoints   = 0
    allocate(gpid(num_gpoints))
endif

allocate(phia_gr_indiv(numnp,num_gpoints))
!-----------------------------------------------------------------------------------------------------------!
end subroutine init_delta
