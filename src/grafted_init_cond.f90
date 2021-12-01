subroutine grafted_init_cond(ns, numnp, gp_filename, qgr, qgr_final)
!------------------------------------------------------------------------------------------------------!
use error_handing
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer, intent(in) :: ns, numnp
integer             :: i1, num_gpoints, gnode_id, iog

character(20), intent(in) :: gp_filename

real(8), intent(out), dimension(numnp,2)    :: qgr
real(8), intent(out), dimension(numnp,ns+1) :: qgr_final
real(8)                                     :: initValue = 0.d0
!------------------------------------------------------------------------------------------------------!
iog = 19
INQUIRE(FILE=gp_filename, EXIST=FILE_EXISTS)

if (FILE_EXISTS) then
    open(unit=iog, file = gp_filename)
else
    write(ERROR_MESSAGE,'(''File '',A15,'' does not exist!'')') gp_filename
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

read(iog,*)
read(iog,*)
read(iog,*)
read(iog,'(I10)') num_gpoints
read(iog,*)
read(iog,*)
read(iog,*)
read(iog,*)
read(iog,*)

!specify grafting points
do i1 = 1, num_gpoints

    read(iog,*) gnode_id, initValue !,x_graft, y_graft, z_graft

    qgr(gnode_id,1)       = initValue
    qgr_final(gnode_id,1) = initValue

    if (gnode_id > numnp) then
        write(ERROR_MESSAGE,'(''ID of grafted chain ('',I10,'') is larger from numnp ('',I10,'')'')') gnode_id, numnp
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif
enddo

close(iog)
!------------------------------------------------------------------------------------------------------!
end subroutine grafted_init_cond
