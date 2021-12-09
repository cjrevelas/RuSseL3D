subroutine import_mesh()
!--------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use error_handing
use write_helper, only: adjl
use parser_vars,  only: iow, n_dirichlet_faces, n_nanopart_faces, prof_dim,     &
&                       ids_dirichlet_faces, ids_nanopart_faces
use geometry,     only: ix, xc, numnp, numel, el_node, is_dir_face, max_el_node,&
&                       node_in_q0_face, nel, n_el_node, ndm, con_l2, all_el,   &
&                       box_lo, box_hi, box_len
use kcw,          only: F_m
use iofiles,      only: mesh_filename, dir_faces, com_12, inter, mesh_out, mesh_prof
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
character(len=200) :: line, line_internal

integer                              :: nbin, ibin
integer                              :: reason, reason_internal, idummy, i, j, i1, j1, k1, m1, i_node, ii, jj
integer                              :: num_nodes_per_face_elem, num_face_elem
integer, allocatable, dimension(:)   :: face_entity_id, temp3
integer, allocatable, dimension(:,:) :: face_node_id
integer, allocatable, dimension(:)   :: prof_1D_node

real(8) :: box_volume = 0.d0, tol = 1.e-8, prof_bin

!fhash module variables
type(fhash_type__ints_double) :: h
type(ints_type)               :: key
integer                       :: n_keys
integer                       :: key_value
logical                       :: success
!--------------------------------------------------------------------!
inquire(file=mesh_filename, exist=file_exists)

if (file_exists) then
    open(unit=12, file=mesh_filename)
else
    write(ERROR_MESSAGE,'("File ",A15," does not exist!")') mesh_filename
    call exit_with_error(1,1,1,ERROR_MESSAGE)
endif

do
    read(12,'(A100)',IOSTAT=reason) line

    if (reason>0) then
        write(*,*) "Something went wrong!"
    elseif (reason<0) then
        exit
    else
        if (INDEX(line,"# sdim") > 0) then
            read(line,*) ndm
        elseif (INDEX(line,"# number of mesh points") > 0) then
            read(line,*) numnp
            allocate(xc(ndm,numnp))
        elseif (INDEX(line,"# Mesh point coordinates") > 0) then

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

            write(iow,*)
            write(*,*)
            write(iow,'(A85)')adjl('-------------------------------------BOX DIMENSIONS----------------------------------',85)
            write(*  ,'(A85)')adjl('-------------------------------------BOX DIMENSIONS----------------------------------',85)

            write(iow,'(A6,A13,A17,A18)') "dim", "length", "min", "max"
            write(6  ,'(A6,A13,A17,A18)') "dim", "length", "min", "max"
            do j = 1, ndm
                box_len(j) = box_hi(j) - box_lo(j)
                write(iow,'(I5,2X,3(E16.9,2X))')j, box_len(j), box_lo(j), box_hi(j)
                write(6  ,'(I5,2X,3(E16.9,2X))')j, box_len(j), box_lo(j), box_hi(j)
            enddo

            write(6,*)
            write(iow,*)

            box_volume = box_len(1) * box_len(2) * box_len(3)
            write(iow,'(3X,A40,E16.9,A13)')adjl("Box volume:",40),box_volume," [Angstrom^3]"
            write(6  ,'(3X,A40,E16.9,A13)')adjl("Box volume:",40),box_volume," [Angstrom^3]"
        elseif (INDEX(line," tet") > 0) then

            read(12,*)
            read(12,*)
            read(12,*) nel
            read(12,*) numel
            allocate(ix(nel,numel))

            all_el = nel * nel * numel

            write(iow,*)
            write(*,*)
            write(iow,'(A85)')adjl('--------------------------------------SPATIAL MESH-----------------------------------',85)
            write(*  ,'(A85)')adjl('--------------------------------------SPATIAL MESH-----------------------------------',85)

            write(iow,'(3X,"Number of mesh points (numnp):         ",I16)') numnp
            write(iow,'(3X,"Number of elements (numel):            ",I16)') numel
            write(iow,'(3X,"Number of nodes per element (nel):     ",I16)') nel
            write(iow,'(3X,"Number of matrix indeces:              ",I16)') all_el

            write(6,'(3X,"Number of mesh points (numnp):         ",I16)') numnp
            write(6,'(3X,"Number of elements (numel):            ",I16)') numel
            write(6,'(3X,"Number of nodes per element (nel):     ",I16)') nel
            write(6,'(3X,"Number of matrix indeces:              ",I16)') all_el

            read(12,*)

            do i = 1, numel
                read(12,*) (ix(j,i), j = 1, nel)
            enddo

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
        elseif (INDEX(line," tri") > 0) then

            read(12,*)
            read(12,*)
            read(12,*) num_nodes_per_face_elem
            read(12,*) num_face_elem
            allocate(face_node_id(num_nodes_per_face_elem,num_face_elem))

            read(12,*)

            do i = 1, num_face_elem
                read(12,*) (face_node_id(j,i), j = 1, num_nodes_per_face_elem)
            enddo

            face_node_id = face_node_id + 1

            allocate(face_entity_id(num_face_elem))

            do
                read(12,'(A100)',IOSTAT=reason_internal) line_internal

                if (reason_internal>0) then
                    write(*,*) "Something went wrong!"
                elseif (reason_internal<0) then
                    write(*,*) "Yaouza"
                    exit
                else
                    if(INDEX(line_internal,'# number of geometric entity indices') > 0) then
                        read(12,*)

                        do i = 1, num_face_elem
                            read(12,*)  face_entity_id(i)
                        enddo

                        exit
                    endif
                endif
            enddo

            face_entity_id = face_entity_id + 1
        endif
    endif
enddo

close(12)

!compute the maximum number of elements per node

! compute the array elements of node

allocate(n_el_node(numnp))
n_el_node = 0

max_el_node = 0
do ii = 1, numel
   do jj = 1, 4
      i_node = ix(jj, ii)
      n_el_node(i_node) = n_el_node(i_node) + 1
      max_el_node = MAX(max_el_node, n_el_node(i_node))
!      write(*,*)max_el_node, n_el_node(i_node)
   enddo
enddo

write(iow,'(3X,"Maximum number of elements per node: ",I18)') max_el_node
write(6  ,'(3X,"Maximum number of elements per node: ",I18)') max_el_node

allocate(el_node(numnp, max_el_node))
el_node = 0

n_el_node = 0
do ii = 1, numel
   do jj = 1, 4
      i_node = ix(jj, ii)
      n_el_node(i_node) = n_el_node(i_node) + 1
!      write(*,*)max_el_node, n_el_node(i_node)
      el_node(i_node, n_el_node(i_node)) = ii
   enddo
enddo

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
allocate(key%ints(2)) !each key is defined by a pair (2) of nodes

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
        enddo
    enddo
enddo
call h%clear()

!determine all elements belonging to dirichlet faces
allocate(node_in_q0_face(numnp))
node_in_q0_face = .false.

#ifdef DEBUG_OUTPUTS
open(unit=123, file = dir_faces)
#endif

is_dir_face = .false.

do j = 1, num_face_elem
    do i1 = 1, n_dirichlet_faces
        if (face_entity_id(j)==ids_dirichlet_faces(i1)) then
            do i = 1, num_nodes_per_face_elem
                idummy = face_node_id(i,j)

                !if the node belongs to a dirichlet face
                node_in_q0_face(idummy) = .True.

                !find if a node is located at a corner
                k1 = 0
                do m1 = 1, ndm
                    if (DABS(xc(m1, idummy) - box_lo(m1)) < tol) then
                        k1 = k1 + 1
                    endif
                    if (DABS(xc(m1, idummy) - box_hi(m1)) < tol) then
                        k1 = k1 + 1
                    endif
                enddo

                !if a node is located at a corner skip the loop
                if (k1 > 1) cycle

                do m1 = 1, ndm
                    if (DABS(xc(m1, idummy) - box_lo(m1)) < tol) then
                        is_dir_face(m1,1) = .true.
                    endif
                    if (DABS(xc(m1, idummy) - box_hi(m1)) < tol) then
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
                node_in_q0_face(idummy) = .True.
            enddo
        endif
    enddo
enddo
!-------------------------DEBUG SECTION-----------------------------!
#ifdef DEBUG_OUTPUTS
close(123)
close(13)

open(unit=77, file = com_12)
do i1 = 1, all_el
     write(77,'(4I8,L8)') i1, F_m%row(i1), F_m%col(i1), con_l2(i1), con_l2(i1)/=i1
enddo
close(77)

open (unit=77, file = inter)
do i = 1, numel
    write(77,'(I9,10(2X,I9))') (ix(j,i), j = 1, nel)
enddo
close(77)

open(77, file = mesh_out)
write(77,'(A13,2A22)') 'x', 'y', 'z'
do i = 1, numnp
    write(77,'(F21.15,1X,F21.15,1X,F21.15)') (xc(j,i), j = 1, ndm)
enddo
close(77)
!-------------------------------------------------------------------!
! APS profile section
! APS TEMP: this should be encapsulated into a subroutine
prof_bin = 0.5d0
nbin = NINT((box_hi(prof_dim) - box_lo(prof_dim)) / prof_bin) + 1

allocate(prof_1D_node(nbin))

prof_1D_node=0
do i = 1, numnp
    ibin = NINT((xc(prof_dim,i) - box_lo(prof_dim))/prof_bin) + 1
    prof_1D_node(ibin) = prof_1D_node(ibin) + 1
enddo

open(77, file = mesh_prof)
do i = 1, nbin
    write(77,*)i, prof_1D_node(i)
enddo
close(77)
#endif

return
!--------------------------------------------------------------------!
end subroutine import_mesh
