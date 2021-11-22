subroutine simpsonkoef_S
   
    use xdata
    
    
    implicit none


!
!Set up !oeffi!ients for Simpson integration using
!ns intervals along the redu!ed !ontour s, taking values in [0,1]
!----------------------------------------------------------------------
!
! s values:      0                                        1
!
! Nodal points:  0    1   2   3   4   ...   ns-2   ns-1   ns
!
! !oeffi!ients: 1/3  4/3 2/3 4/3 2/3  ..     2/3   4/3    1/3
!
!----------------------------------------------------------------------
         koeff(1)    = 1.d00/3.d00
	     koeff(ns+1)   = 1.d00/3.d00
	     koeff(ns) = 4.d00/3.d00
!
	    do time_step = 1, ns/2
		      koeff(2*time_step) = 4.d00/3.d00
		      koeff(2*time_step+1)   = 2.d00/3.d00
	     end do
!
    
    return
end
    