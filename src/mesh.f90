subroutine mesh
!--------------------------------------------------------------------!
!fhash module
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
!/fhash module
use error_handing
use write_helper
use parser_vars
use kcw
use geometry
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
character(len=15) :: dummy
character(len=14) :: mesh_filename

integer :: idummy, i, j, i1, j1, k1, m1
integer :: num_elem_types

integer :: num_nodes_per_vertex_elem, num_vertex_elem
integer :: num_params_per_vertex_elem, num_vertex_params
integer :: num_nodes_per_edge_elem, num_edge_elem
integer :: num_params_per_edge_elem, num_edge_params
integer :: num_nodes_per_face_elem, num_face_elem
integer :: num_params_per_face_elem, num_face_params

integer, allocatable, dimension(:)   :: vertex_entity_id, edge_entity_id, face_entity_id, temp3
integer, allocatable, dimension(:,:) :: vertex_node_id, edge_node_id, face_node_id
real(8), allocatable, dimension(:,:) :: vertex_param, edge_param, face_param

real(8) :: tol = 1.e-8

!profile section variables
integer, allocatable, dimension(:) :: prof_1D_node
real(8)                            :: prof_bin
integer                            :: nbin, ibin

!fhash module variables
type(fhash_type__ints_double) :: h
type(ints_type) :: key
integer :: n_keys
integer :: key_value
logical :: success
!--------------------------------------------------------------------!
mesh_filename = 'mesh.in.mphtxt'

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

!reading meshpoint coordinates and calculating box dimensions
box_lo  = 0.d0
box_hi  = 0.d0
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

!read types and numbers of elements
read (12,'(A60)') dummy

read (12,*) num_elem_types

do i = 1, 6
    read (12,'(A60)') dummy
enddo

!----------------vertex elements------------------!
read(12,*) num_nodes_per_vertex_elem
read(12,*) num_vertex_elem
allocate(vertex_node_id(num_nodes_per_vertex_elem,num_vertex_elem))

call read_integer_value(num_nodes_per_vertex_elem, num_vertex_elem, vertex_node_id)

read (12,'(A60)') dummy

read(12,*) num_params_per_vertex_elem
read(12,*) num_vertex_params
allocate(vertex_param(num_params_per_vertex_elem,num_vertex_params))

call read_real_value(num_params_per_vertex_elem, num_vertex_params, vertex_param)

read (12,'(A60)') dummy

allocate(vertex_entity_id(num_vertex_elem))

call entity(num_vertex_elem, vertex_entity_id)

do i = 1, 8
    read (12,'(A60)') dummy
enddo

!------------------edge elements------------------!
read(12,*) num_nodes_per_edge_elem
read(12,*) num_edge_elem
allocate(edge_node_id(num_nodes_per_edge_elem,num_edge_elem))

call read_integer_value(num_nodes_per_edge_elem, num_edge_elem, edge_node_id)

read(12,'(A60)') dummy

read(12,*) num_params_per_edge_elem
read(12,*) num_edge_params
allocate(edge_param(num_params_per_edge_elem,num_edge_params))

call read_real_value(num_params_per_edge_elem, num_edge_params, edge_param)

read (12,'(A60)') dummy

allocate(edge_entity_id(num_edge_elem))

call entity(num_edge_elem, edge_entity_id)

do i = 1, 9
    read (12,'(A60)') dummy
enddo

!------------------face elements------------------!
read(12,*) num_nodes_per_face_elem
read(12,*) num_face_elem
allocate(face_node_id(num_nodes_per_face_elem,num_face_elem))

call read_integer_value(num_nodes_per_face_elem, num_face_elem, face_node_id)

read (12,'(A60)') dummy

face_node_id = face_node_id + 1 !change the values of face elements so that they start from 1 instead from 0

read(12,*) num_params_per_face_elem
read(12,*) num_face_params
allocate(face_param(num_params_per_face_elem,num_face_params))

call read_real_value(num_params_per_face_elem, num_face_params, face_param)

read (12,'(A60)') dummy

allocate(face_entity_id(num_face_elem))

call entity(num_face_elem, face_entity_id)

face_entity_id = face_entity_id + 1 !change the value of face id's so that they start from 1 instead from 0

!-------------------------------------------------!
!up/down pairs
read(12,'(A60)') dummy
read(12,*) idummy

do i = 1, idummy
    read(12,'(A60)') dummy
enddo

do i = 1, 6
    read(12,'(A60)') dummy
enddo

!------------------domain elements------------------!
read(12,*) nel
read(12,*) numel
allocate(ix(nel,numel))

call read_integer_value(nel, numel, ix)

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
!---------------------------------------------------!
!allocate and initialize arrays for matrix assembly
allocate(F_m%row(all_el))
allocate(F_m%col(all_el))
allocate(F_m%g(all_el))
allocate(F_m%rh(all_el))
allocate(F_m%c(all_el))
allocate(F_m%k(all_el))
allocate(F_m%w(all_el))

F_m%g   = 0.d0
F_m%k   = 0.d0
F_m%c   = 0.d0
F_m%rh  = 0.d0
F_m%row = 0
F_m%col = 0

allocate(con_l2(all_el))
con_l2 = 0

!assembly the con_12 (hash) matrix
allocate(key%ints(2)) !each key points to a pair (2) of nodes

!total number of required keys
n_keys = nel * numel
call h%reserve(n_keys*2)

all_el = 0
do m1 = 1, numel
    do j1 = 1, nel
        do i1 = 1, nel
            all_el = all_el + 1

            F_m%row(all_el) = ix(j1,m1)
            F_m%col(all_el) = ix(i1,m1)

            !define the pair of nodes to be examined and assigned a key_value
            key%ints(1) = ix(j1,m1)
            key%ints(2) = ix(i1,m1)

            !assign key_value to the pair
            call h%get(key, key_value, success)
            if (success) then
               con_l2(all_el) = key_value      !this pair has already been met, thus assigned a key_value
            else
               call h%set(key, all_el)         !store the new key_value for next iteration's check
               con_l2(all_el) = all_el         !this pair is met for the first time
            endif
        end do
    end do
end do
call h%clear()

!determine all elements belonging to dirichlet faces
allocate(elem_in_q0_face(numnp))
elem_in_q0_face = .false.

write(iow,'(/A54)')'* Find all elements belonging to dirichlet q=0 faces..'
write(6  ,'(/A54)')'* Find all elements belonging to dirichlet q=0 faces..'

#ifdef DEBUG_OUTPUTS
open(unit=123, file='dir_faces.out.txt')
#endif

is_dir_face      = .false.

do j = 1, num_face_elem
    do i1 = 1, n_dirichlet_faces
        if (face_entity_id(j)==ids_dirichlet_faces(i1)) then
            do i = 1, num_nodes_per_face_elem
                idummy= face_node_id(i,j)

                !if the node belongs to a dirichlet face
                elem_in_q0_face(idummy) = .True.

                !find if a node is located at a corner
                k1 = 0
                do m1 = 1, ndm
                    if (dabs(xc(m1, idummy) - box_lo(m1)) < tol) then
                        k1 = k1 + 1
                    endif
                    if (dabs(xc(m1, idummy) - box_hi(m1)) < tol) then
                        k1 = k1 + 1
                    endif
                enddo

                !if a node is located at a corner skip the loop
                if (k1 > 1) cycle

                do m1 = 1, ndm
                    if (dabs(xc(m1, idummy) - box_lo(m1)) < tol) then
                        is_dir_face(m1,1) = .true.
                    endif
                    if (dabs(xc(m1, idummy) - box_hi(m1)) < tol) then
                        is_dir_face(m1,2) = .true.
                    endif
                enddo
#ifdef DEBUG_OUTPUTS
                write(123,'(3(E16.8),6(L3))') xc(1, idummy), xc(2, idummy), xc(3, idummy), &
    &                                 (is_dir_face(k1,1), k1 = 1,3), (is_dir_face(k1,2), k1 = 1,3)
#endif
            enddo
        endif
    enddo

    !nanoparticles section
    do i1 = 1, n_nanopart_faces
        if (face_entity_id(j)==ids_nanopart_faces(i1)) then
            do i = 1, num_nodes_per_face_elem
                idummy = face_node_id(i,j)

                !if the node belongs to a dirichlet face
                elem_in_q0_face(idummy) = .True.
            enddo
        endif
    enddo
enddo
!-------------------------DEBUG SECTION-----------------------------!
#ifdef DEBUG_OUTPUTS
close(123)
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
!-------------------------------------------------------------------!
! APS profile section
! APS TEMP: this should be encapsulated into a subroutine
prof_bin = 0.5d0
nbin = nint((box_hi(prof_dim) - box_lo(prof_dim)) / prof_bin) + 1

allocate(prof_1D_node(nbin))

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

return
!--------------------------------------------------------------------!
end subroutine mesh



subroutine read_integer_value(num_per_elem, num_elem, id)
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer :: num_per_elem, num_elem, i, j

integer, dimension(num_per_elem, num_elem) :: id

character(len=60) :: dum
!--------------------------------------------------------------------!
read (12,'(A60)') dum

if (num_elem.ne.0) then
    if (num_per_elem.eq.1) then
        do i = 1, num_elem
            read(12,*)  id(1,i)
        enddo
    else
        do i = 1, num_elem
            read(12,*)  (id(j,i), j = 1, num_per_elem)
        enddo
    endif
endif

return
!--------------------------------------------------------------------!
end subroutine read_integer_value



subroutine entity(num_elem, entity_id)
!--------------------------------------------------------------------!
use error_handing
implicit none
!--------------------------------------------------------------------!
integer :: num_elem, num_entities, i

integer, dimension(num_elem) :: entity_id

character(len=60) :: dum
!--------------------------------------------------------------------!
read(12,*) num_entities

read (12,'(A60)') dum

if (num_entities.eq.num_elem) then
    do i = 1, num_entities
        read(12,*)  entity_id(i)
    enddo
else
    ERROR_MESSAGE='Error in entity..'
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

return
!--------------------------------------------------------------------!
end subroutine entity



subroutine read_real_value(num_per_elem, num_params, id)
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
integer:: num_per_elem, num_params, i, j

real(8), dimension(num_per_elem, num_params) :: id

character(len=60) :: dum
!--------------------------------------------------------------------!
read (12,'(A60)') dum

if (num_params.ne.0) then
    if (num_per_elem.eq.1) then
        do i = 1, num_params
            read(12,*)  id(1,i)
        enddo
    else
        do i = 1, num_params
            read(12,*)  (id(j,i), j = 1, num_per_elem)
        enddo
    endif
endif

return
!--------------------------------------------------------------------!
end subroutine read_real_value
