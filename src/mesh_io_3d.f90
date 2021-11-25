subroutine mesh_io_3d
!--------------------------------------------------------------------!
! APS 17/08/19: ADD. fhash module
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
!/ fhash module
use mdata
use xdata
use kcw
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!    
character(len=15) :: dummy

integer :: idummy
integer :: i, sdim, j
integer :: nmeltypes
integer :: vtxnum, vtxel
integer :: vtxparnum, vtxparel
integer :: edgnum, edgel
integer :: edgparnum, edgparel
integer :: fcparnum, fcparel

integer, allocatable, dimension(:,:) :: vtxelement, edgelement
integer, allocatable, dimension(:,:) :: vtxpar 
integer, allocatable, dimension(:) :: vtxentity, edgentity, temp3

real(8), allocatable, dimension(:,:) :: edgpar, fcpar

! APS 17/08/19: ADD. fhash module
type(fhash_type__ints_double) :: h
type(ints_type) :: key
integer :: n_keys
integer :: key_value
logical :: success

!/ fhash module
!--------------------------------------------------------------------!
open(unit=12, file='Fem_3D.mphtxt')

do i = 1, 17
    read(12,'(A60)') dummy
enddo

read(12,*) sdim
read(12,*) numnp


allocate(xc(sdim,numnp))

do i = 1, 3
    read(12,'(A60)') dummy
enddo

box_lo = 0.d0
box_hi = 0.d0
box_len = 0.d0

do i = 1, numnp
    read(12,*) (xc(j,i), j = 1, sdim)

    do j = 1, sdim
       box_hi(j) = max(xc(j,i), box_hi(j))
       box_lo(j) = min(xc(j,i), box_lo(j))
    enddo
enddo

do j = 1, sdim
    box_len(j) = box_hi(j) - box_lo(j)
    write(6,'(A15,I3,3F20.7)')"box Len lo / hi",j, box_len(j), box_lo(j), box_hi(j)
enddo

volume = box_len(1)*box_len(2)*box_len(3)
write(6,*)"Volume = ", volume

read (12,'(A60)') dummy

!element types
read (12,*) nmeltypes


#ifdef DEBUG_OUTPUTS
open(unit=13, file='com.out.txt')
write(13,'(I9)') sdim
write(13,'(I9)') numnp
do i = 1, numnp
    write(13,'(F19.15,1X,F19.15,1X,F19.15)') (xc(j,i), j = 1, sdim)
enddo
write(13,*) nmeltypes
#endif


do i = 1, 6
    read (12,'(A60)') dummy
enddo

call allocate_matrix(vtxnum, vtxel)

allocate(vtxelement(vtxnum,vtxel))

call read_integer(vtxnum, vtxel, vtxelement)

read (12,'(A60)') dummy

call allocate_matrix(vtxparnum, vtxparel)

allocate(vtxpar(vtxparnum,vtxparel))

call read_integer(vtxparnum, vtxparel, vtxpar)

read (12,'(A60)') dummy

allocate(vtxentity(vtxel))

call entity(vtxel, vtxentity)

do i = 1, 8
    read (12,'(A60)') dummy
enddo

!EDGE ELEMENTS
call allocate_matrix(edgnum, edgel)

allocate(edgelement(edgnum,edgel))

call read_integer(edgnum, edgel, edgelement)

read(12,'(A60)') dummy

call allocate_matrix(edgparnum, edgparel)

allocate(edgpar(edgparnum,edgparel))

call read_real(edgparnum, edgparel, edgpar)

read (12,'(A60)') dummy

allocate(edgentity(edgel))

call entity(edgel, edgentity)

do i = 1, 9
    read (12,'(A60)') dummy
enddo

!FACE ELEMENTS
call allocate_matrix(fcnum, fcel)

allocate(fcelement(fcnum,fcel))

call read_integer(fcnum, fcel, fcelement) 

read (12,'(A60)') dummy

!changes the values of face element so as to start from 1 instead from 0
fcelement = fcelement + 1

call allocate_matrix(fcparnum, fcparel)

allocate(fcpar(fcparnum,fcparel))

call read_real(fcparnum, fcparel, fcpar)

read (12,'(A60)') dummy

allocate(fcentity(fcel))

call entity(fcel, fcentity)
! changes the values of  element entinty so as to start from 1 instead from 0
fcentity = fcentity + 1

!UP/DOWN PAIRS
read(12,'(A60)') dummy
read(12,*) idummy

do i = 1, idummy
    read(12,'(A60)') dummy
enddo

do i = 1, 6
    read(12,'(A60)') dummy
enddo

!DOMAIN ELEMENTS
call allocate_matrix(dmnum, dmel)

allocate(ix(dmnum,dmel))

call read_integer(dmnum, dmel, ix)

do i = 1, dmel
    do j = 1, dmnum
        ix(j,i) = ix(j,i) + 1
    enddo
enddo

ndm   = sdim
nel   = dmnum
numel = dmel

if (nel>4) then
    allocate(temp3(numel))
    do i = 1, numel
        temp3(i) = ix(7,i)
        ix(7,i)  = ix(6,i)
        ix(6,i)  = temp3(i)
    enddo
endif

all_el = nel * nel * numel

allocate(g_m%row(all_el))
allocate(g_m%col(all_el))
allocate(g_m%value(all_el))

allocate(rh_m%row(all_el))
allocate(rh_m%col(all_el))
allocate(rh_m%value(all_el))

allocate(c_m%row(all_el))
allocate(c_m%col(all_el))
allocate(c_m%value(all_el))

allocate(k_m%row(all_el))
allocate(k_m%col(all_el))
allocate(k_m%value(all_el))

allocate(w_m%row(all_el))
allocate(w_m%col(all_el))
allocate(w_m%value(all_el))

allocate(con_l2(all_el))

g_m%value = 0.



! APS 17/08/19: ADD. fhash module

allocate(key%ints(2))
n_keys = nel * numel
call h%reserve(n_keys*2)

all_el = 0
do m1 = 1, numel
    do j1 = 1, nel
        do i1 = 1, nel
            all_el = all_el + 1

            g_m%row(all_el) = ix(j1,m1)
            g_m%col(all_el) = ix(i1,m1)

            key%ints(1) =  ix(j1,m1)
            key%ints(2) =  ix(i1,m1)
            call h%get(key, key_value, success)
            if (success) then
               con_l2(all_el) = key_value
            else
               call h%set(key, all_el)
               con_l2(all_el) = all_el
            endif
        end do
    end do
end do
call h%clear()

!/ new section

! old section
!all_el = 0
!count  = 1
!do m1 = 1, numel
!    do j1 = 1, nel
!        do i1 = 1, nel
!            all_el = all_el + 1
!            g_m%row(all_el) = ix(j1,m1)
!            g_m%col(all_el) = ix(i1,m1)
!            do k1 = 1, all_el-1
!                if ((g_m%row(k1)==g_m%row(all_el)).and.(g_m%col(k1)==g_m%col(all_el))) then
!                    con_l2(all_el) = k1
!                    exit
!                elseif (k1==all_el-1) then 
!                    count = count + 1
!                    con_l2(all_el) = all_el
!               endif
!            enddo
!        end do
!    end do
!end do
!con_l2(1)=1
!/ old section

c_m%row  = g_m%row
k_m%row  = g_m%row
w_m%row  = g_m%row
rh_m%row = g_m%row
c_m%col  = g_m%col
k_m%col  = g_m%col
w_m%col  = g_m%col
rh_m%col = g_m%col

g_m%value = 0.
k_m%value = 0.
c_m%value = 0.

allocate(rdiag1(numnp))

#ifdef DEBUG_OUTPUTS
close(13)

open(unit=77, file='com_12.out.txt')
do i1 = 1, all_el
     write(77,*) i1, g_m%row(i1), g_m%col(i1), con_l2(i1), con_l2(i1)/=i1
enddo
close(77)

open (unit=77, file='inter.out.txt')  
do i = 1, dmel
    write(77,'(I9,10(2X,I9))') (ix(j,i), j = 1, dmnum) 
enddo
close(77)

open(77, file = 'mesh.out.txt')
do i = 1, numnp
    write(77,'(F21.15,1X,F21.15,1X,F21.15)') (xc(j,i), j = 1, sdim) 
enddo
close(77)
#endif

return
!--------------------------------------------------------------------!
end subroutine  mesh_io_3d



subroutine allocate_matrix(num, el)
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: num, el
!--------------------------------------------------------------------!
read(12,*) num
read(12,*) el

return
!--------------------------------------------------------------------!
end subroutine allocate_matrix



subroutine read_integer(num, el, element)
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: num, el, i, j

integer :: element(num,el)

character(len=60) :: dum
!--------------------------------------------------------------------!
read (12,'(A60)') dum

if (el.ne.0) then
    if (num.eq.1) then
        do i = 1, el
            read(12,*)  element(1,i)
#ifdef DEBUG_OUTPUTS
            write(13,*) element(1,i)
#endif
        enddo
    else
        do i = 1, el
            read(12,*)  (element(j,i), j = 1, num)
#ifdef DEBUG_OUTPUTS
            write(13,*) (element(j,i), j = 1, num)
#endif
        enddo
    endif
endif

return
!--------------------------------------------------------------------!
end subroutine read_integer

subroutine entity(el, endity)
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: el, equal, i!, num

integer :: endity(el)

character(len=60) :: dum
!--------------------------------------------------------------------!
read(12,*)equal

read (12,'(A60)') dum

if (equal.eq.el) then
    do i = 1, el
        read(12,*)  endity(i)
#ifdef DEBUG_OUTPUTS
        write(13,*) endity(i)
#endif
    enddo
else
    write(6,*) "Error in entity"
endif

return
!--------------------------------------------------------------------!
end subroutine entity



subroutine read_real(num, el, element)
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer:: num, el, i, j

real(8) :: element(num,el)

character(len=60) :: dum
!--------------------------------------------------------------------!
read (12,'(A60)') dum

if (el.ne.0) then
    if (num.eq.1) then
        do i = 1, el
            read(12,*)  element(1,i)
#ifdef DEBUG_OUTPUTS
            write(13,*) element(1,i)
#endif
        enddo
    else
        do i = 1, el
            read(12,*)  (element(j,i), j = 1, num)
#ifdef DEBUG_OUTPUTS
            write(13,*) (element(j,i), j = 1, num)
#endif
        enddo
    endif
endif

return
!--------------------------------------------------------------------!
end subroutine read_real
