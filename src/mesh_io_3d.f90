subroutine mesh_io_3d
!--------------------------------------------------------------------!
!fhash module
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
!/fhash module
use error_handing
use write_helper
use xdata
use kcw
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
character(len=15) :: dummy
character(len=13) :: mesh_filename

integer                              :: idummy
integer                              :: i, j, i1, j1, k1, m1
integer                              :: nmeltypes
integer                              :: vtxnum, vtxel
integer                              :: vtxparnum, vtxparel
integer                              :: edgnum, edgel
integer                              :: edgparnum, edgparel
integer                              :: fcparnum, fcparel
integer                              :: fcnum, fcel
integer, allocatable, dimension(:,:) ::  fcelement
integer, allocatable, dimension(:)   :: fcentity
integer, allocatable, dimension(:,:) :: vtxelement, edgelement
integer, allocatable, dimension(:,:) :: vtxpar
integer, allocatable, dimension(:)   :: vtxentity, edgentity, temp3

real(8), allocatable, dimension(:,:) :: edgpar, fcpar

!profile section
integer, allocatable, dimension(:) :: prof_1D_node
real(8)                            :: prof_bin
integer                            :: nbin, ibin

!fhash module
type(fhash_type__ints_double) :: h
type(ints_type) :: key
integer :: n_keys
integer :: key_value
logical :: success
!/ fhash module
!--------------------------------------------------------------------!
mesh_filename = 'Fem_3D.mphtxt'

write(iow,'(/A42,A15)')adjl('*Reading mesh from file:',42),mesh_filename
write(6  ,'(/A42,A15)')adjl('*Reading mesh from file:',42),mesh_filename

INQUIRE(FILE=mesh_filename, EXIST=FILE_EXISTS)

if (FILE_EXISTS) then
    open(unit=12, file=mesh_filename)
else
    write(ERROR_MESSAGE,'(''File '',A15,'' does not exist!'')')mesh_filename
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

do i = 1, 17
    read(12,'(A60)') dummy
enddo

read(12,*) ndm
read(12,*) numnp


allocate(xc(ndm,numnp))

do i = 1, 3
    read(12,'(A60)') dummy
enddo

!*******************************************************************!
!                Reading Mesh points and box dimensions             !
!*******************************************************************!
box_lo = 0.d0
box_hi = 0.d0
box_len = 0.d0

do i = 1, numnp
    read(12,*) (xc(j,i), j = 1, ndm)

    do j = 1, ndm
       box_hi(j) = max(xc(j,i), box_hi(j))
       box_lo(j) = min(xc(j,i), box_lo(j))
    enddo
enddo

write(iow,'(/''Box dimensions..'')')
write(iow,'(A5,2x,3(A16,2x))') 'dim','length','min','max'
write(6  ,'(/''Box dimensions..'')')
write(6  ,'(A5,2x,3(A16,2x))') 'dim','length','min','max'
do j = 1, ndm
    box_len(j) = box_hi(j) - box_lo(j)
    write(iow,'(I5,2x,3(E16.9,2x))')j, box_len(j), box_lo(j), box_hi(j)
    write(6  ,'(I5,2x,3(E16.9,2x))')j, box_len(j), box_lo(j), box_hi(j)
enddo

volume = box_len(1)*box_len(2)*box_len(3)
write(iow,'(/A43,E16.9,A11)')adjl('System volume:',43),volume,' Angstrom^3'
write(6  ,'(/A43,E16.9,A11)')adjl('System volume:',43),volume,' Angstrom^3'

!*******************************************************************!
!                       Deal with the elements                      !
!*******************************************************************!
read (12,'(A60)') dummy

!element types
read (12,*) nmeltypes

do i = 1, 6
    read (12,'(A60)') dummy
enddo

!VERTEX ELEMENTS
call read_matrix_bounds(vtxnum, vtxel)

allocate(vtxelement(vtxnum,vtxel))

call read_integer(vtxnum, vtxel, vtxelement)

read (12,'(A60)') dummy

call read_matrix_bounds(vtxparnum, vtxparel)

allocate(vtxpar(vtxparnum,vtxparel))

call read_integer(vtxparnum, vtxparel, vtxpar)

read (12,'(A60)') dummy

allocate(vtxentity(vtxel))

call entity(vtxel, vtxentity)

do i = 1, 8
    read (12,'(A60)') dummy
enddo

!EDGE ELEMENTS
call read_matrix_bounds(edgnum, edgel)

allocate(edgelement(edgnum,edgel))

call read_integer(edgnum, edgel, edgelement)

read(12,'(A60)') dummy

call read_matrix_bounds(edgparnum, edgparel)

allocate(edgpar(edgparnum,edgparel))

call read_real(edgparnum, edgparel, edgpar)

read (12,'(A60)') dummy

allocate(edgentity(edgel))

call entity(edgel, edgentity)

do i = 1, 9
    read (12,'(A60)') dummy
enddo

!FACE ELEMENTS
call read_matrix_bounds(fcnum, fcel)

allocate(fcelement(fcnum,fcel))

call read_integer(fcnum, fcel, fcelement)

read (12,'(A60)') dummy

!change the values of face element so as to start from 1 instead from 0
fcelement = fcelement + 1

call read_matrix_bounds(fcparnum, fcparel)

allocate(fcpar(fcparnum,fcparel))

call read_real(fcparnum, fcparel, fcpar)

read (12,'(A60)') dummy

allocate(fcentity(fcel))

call entity(fcel, fcentity)
!changes the values of element entinty so as to start from 1 instead from 0
fcentity = fcentity + 1

!up/down pairs
read(12,'(A60)') dummy
read(12,*) idummy

do i = 1, idummy
    read(12,'(A60)') dummy
enddo

do i = 1, 6
    read(12,'(A60)') dummy
enddo

!DOMAIN ELEMENTS
call read_matrix_bounds(nel, numel)

allocate(ix(nel,numel))

call read_integer(nel, numel, ix)

do i = 1, numel
    do j = 1, nel
        ix(j,i) = ix(j,i) + 1
    enddo
enddo

if (nel>4) then
    allocate(temp3(numel))
    do i = 1, numel
        temp3(i) = ix(7,i)
        ix(7,i)  = ix(6,i)
        ix(6,i)  = temp3(i)
    enddo
endif

all_el = nel * nel * numel

!*******************************************************************!
!                  Allocate and initialize the arrays               !
!*******************************************************************!
allocate(F_m%row(all_el))
allocate(F_m%col(all_el))
allocate(F_m%g(all_el))
allocate(F_m%rh(all_el))
allocate(F_m%c(all_el))
allocate(F_m%k(all_el))
allocate(F_m%w(all_el))

F_m%g = 0.d0
F_m%k = 0.d0
F_m%c = 0.d0
F_m%rh = 0.d0
F_m%row = 0
F_m%col = 0

allocate(con_l2(all_el))
con_l2 = 0


!*******************************************************************!
!                   construct the con_12 matrix                     !
!*******************************************************************!

allocate(key%ints(2))
n_keys = nel * numel
call h%reserve(n_keys*2)

all_el = 0
do m1 = 1, numel
    do j1 = 1, nel
        do i1 = 1, nel
            all_el = all_el + 1

            F_m%row(all_el) = ix(j1,m1)
            F_m%col(all_el) = ix(i1,m1)

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

!*******************************************************************!
!          Find all elements belonging to dirichlet q=0 faces       !
!*******************************************************************!
allocate(elem_in_q0_face(numnp))
elem_in_q0_face = .false.

write(iow,'(/A54)')'* Find all elements belonging to dirichlet q=0 faces..'
write(6  ,'(/A54)')'* Find all elements belonging to dirichlet q=0 faces..'

do j = 1, fcel
    do i1 = 1, n_dirichlet_faces
        if (fcentity(j)==ids_dirichlet_faces(i1))then
            do i = 1, fcnum
                idummy= fcelement(i,j)
                elem_in_q0_face(idummy) = .True.
            enddo
        endif
    enddo
enddo

!*******************************************************************!
!                          DEBUG SECTION                            !
!*******************************************************************!
#ifdef DEBUG_OUTPUTS
close(13)

open(unit=77, file='com_12.out.txt')
do i1 = 1, all_el
     write(77,*) i1, F_m%row(i1), F_m%col(i1), con_l2(i1), con_l2(i1)/=i1
enddo
close(77)

open (unit=77, file='inter.out.txt')
do i = 1, numel
    write(77,'(I9,10(2X,I9))') (ix(j,i), j = 1, nel)
enddo
close(77)

open(77, file = 'mesh.out.txt')
do i = 1, numnp
    write(77,'(F21.15,1X,F21.15,1X,F21.15)') (xc(j,i), j = 1, ndm)
enddo
close(77)


! APS profile section
! APS TEMP: this should be encapsulated into a subroutine
! APS: DO NOT PANIC
prof_bin = 0.5d0
nbin = nint((box_hi(prof_dim) - box_lo(prof_dim)) / prof_bin) + 1

allocate(prof_1D_node(nbin))

!write(6,*)"Profile options"
!write(6,*)prof_dim, prof_bin, nbin

prof_1D_node=0
do i = 1, numnp
    ibin = nint((xc(prof_dim,i) - box_lo(prof_dim))/prof_bin) + 1
    prof_1D_node(ibin) = prof_1D_node(ibin) + 1
enddo

open(77, file = 'prof_mp.out.txt')
do i = 1, nbin
    write(77,*)i, prof_1D_node(i)
enddo
close(77)
#endif

!APS TODO: move this to a better place!!!
print_ev = 3
open (unit=120, file = 'header.out.txt')
do k1 = 1, numnp, print_ev
    write(120,'(E13.5)',advance='no') xc(prof_dim,k1)
enddo
close(120)

return
!--------------------------------------------------------------------!
end subroutine  mesh_io_3d



subroutine read_matrix_bounds(num, el)
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: num, el
!--------------------------------------------------------------------!
read(12,*) num
read(12,*) el

return
!--------------------------------------------------------------------!
end subroutine read_matrix_bounds



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
        enddo
    else
        do i = 1, el
            read(12,*)  (element(j,i), j = 1, num)
        enddo
    endif
endif

return
!--------------------------------------------------------------------!
end subroutine read_integer

subroutine entity(el, endity)
!--------------------------------------------------------------------!
use error_handing
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
    enddo
else
    ERROR_MESSAGE='Error in entity..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
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
        enddo
    else
        do i = 1, el
            read(12,*)  (element(j,i), j = 1, num)
        enddo
    endif
endif

return
!--------------------------------------------------------------------!
end subroutine read_real
