module mdata
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: dmnum, dmel

integer :: numel, nel, ndm, numnp 

integer :: fcnum, fcel

integer, allocatable, dimension(:,:) :: ix, fcelement
integer, allocatable, dimension(:) :: fcentity

real(8), allocatable, dimension(:,:) :: xc
!--------------------------------------------------------------------!	 
end module mdata
