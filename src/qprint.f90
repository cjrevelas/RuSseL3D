subroutine qprint
    use xdata
    use mdata
    
    implicit none
	character(len=20) frmt
	
	
	 WRITE(frmt,'("(3e20.9,",i4,"e20.9)")') ns+1

         open (file='qfree.txt',unit=363)
         do i1=1,numnp
                 !write (363,"(3f20.9,"//ADJUSTR(frmt)//"f20.9)") (xc(ii2,i1),ii2=1,ndm),(qf_final(i1,time_step),time_step=1,ns+1) 
                 write (363,frmt) (xc(ii2,i1),ii2=1,ndm),(qf_final(i1,time_step),time_step=1,ns+1)
         end do  !i1
         close(363)
!

    return 
end
