subroutine init_arrays
!-----------------------------------------------------------------------------------------------------------!
use geometry
use kcw
use parser_vars
use arrays
use delta
use error_handing
use iofiles
!-----------------------------------------------------------------------------------------------------------!
implicit none
!-----------------------------------------------------------------------------------------------------------!
integer :: i1, iog

allocate(rdiag1(numnp))
allocate(wa(numnp),wa_mix(numnp),wa_new(numnp),Ufield(numnp))
allocate(ds_matrix_ed(ns_matrix_ed+1),ds_matrix_conv(ns_matrix_conv+1))
allocate(xs_matrix_ed(ns_matrix_ed+1),xs_matrix_conv(ns_matrix_conv+1))
allocate(koeff_matrix_ed(ns_matrix_ed+1),koeff_matrix_conv(ns_matrix_conv+1))
allocate(phia_mx(numnp))
allocate(qm(numnp,2))
allocate(qm_final(numnp,ns_matrix_ed+1))
allocate(qm_interp_mm(numnp,ns_matrix_conv+1))
allocate(qm_interp_mg(numnp,ns_gr_conv+1))
allocate(qgr_final(numnp,ns_gr_ed+1))
allocate(qgr_interp(numnp,ns_gr_conv+1))
allocate(phia_gr(numnp))

wa                = 0.d0
wa_mix            = 0.d0
wa_new            = 0.d0
Ufield            = 0.d0
rdiag1            = 0.d0
ds_matrix_ed      = 0.d0
ds_matrix_conv    = 0.d0
xs_matrix_ed      = 0.d0
xs_matrix_conv    = 0.d0
koeff_matrix_ed   = 0.d0
koeff_matrix_conv = 0.d0
phia_mx           = 0.d0
qm_final          = 0.d0
qm                = 0.d0
qm_interp_mm      = 0.d0
qm_interp_mg      = 0.d0
qgr_final         = 0.d0
qgr_interp        = 0.d0
phia_gr           = 0.d0

if (use_grafted.eq.1) then
    allocate(ds_gr_ed(ns_gr_ed+1),ds_gr_conv(ns_gr_conv+1))
    allocate(xs_gr_ed(ns_gr_ed+1),xs_gr_conv(ns_gr_conv+1))
    allocate(koeff_gr_ed(ns_gr_ed+1),koeff_gr_conv(ns_gr_conv+1))
    allocate(qgr(numnp,2))

    ds_gr_ed      = 0.d0
    ds_gr_conv    = 0.d0
    xs_gr_ed      = 0.d0
    xs_gr_conv    = 0.d0
    koeff_gr_ed   = 0.d0
    koeff_gr_conv = 0.d0
    qgr           = 0.d0

    iog = 19 

    inquire(file = gp_filename, exist = file_exists)

    if (file_exists) then
        open(unit=iog, file=gp_filename)
    else
        write(ERROR_MESSAGE,'("File ",A15," does not exist!")') gp_filename
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

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

    !specify grafting points
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
!-----------------------------------------------------------------------------------------------------------!
end subroutine init_arrays
