module mdata
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: dmnum, dmel

integer :: numel, nel, ndm, numnp 

integer :: fcnum, fcel
 
integer, allocatable :: ix(:,:), fcelement(:,:), fcentity(:)

real(8), allocatable :: xc(:,:)
!--------------------------------------------------------------------!	 
end module mdata
