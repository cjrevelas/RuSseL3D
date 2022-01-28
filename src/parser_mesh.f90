!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine parser_mesh()
!----------------------------------------------------------------------------------------------------------------------------!
use, intrinsic :: iso_fortran_env
use fhash_module__ints_double
use ints_module
use error_handing_mod
use write_helper_mod, only: adjl
use parser_vars_mod,  only: iow, periodic_axis_id, domain_is_periodic
use geometry_mod,     only: nel, num_of_elems_of_node, box_lo, box_hi, box_len,           &
&                       xc, numnp, numel, global_node_id_type_domain,                     &
&                       ndm, node_pair_id, num_of_bulk_pairs, total_num_of_node_pairs,    &
&                       node_pairing_xx_hash, node_pairing_yy_hash, node_pairing_zz_hash, &
&                       num_dest_xx_neighbors, num_dest_yy_neighbors,                     &
&                       num_dest_zz_neighbors, nen_type_face, numel_type_face
use kcw_mod,          only: F_m
use iofiles_mod,      only: mesh_filename, dir_faces, com_12, inter, mesh_out,         &
&                       mesh_prof, xface1_elements, xface2_elements,                   &
&                       yface1_elements, yface2_elements, zface1_elements,             &
&                       zface2_elements, node_pairing_xx, node_pairing_yy,             &
&                       node_pairing_zz
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
character(len=200) :: line, aux_line

integer                              :: ii, jj, kk, pp, forward_steps
integer                              :: starting_pair, ending_pair
integer                              :: reason, id_of_entity_it_belongs, max_num_of_elems_per_node
integer                              :: nen_type_vertex, numel_type_vertex
integer                              :: nen_type_edge, numel_type_edge
integer, allocatable, dimension(:)   :: temp3
integer, allocatable, dimension(:,:) :: global_node_id_type_vertex, global_node_id_type_edge, global_node_id_type_face

real(8) :: box_volume = 0.d0

logical :: success

type(fhash_type__ints_double) :: elemcon

type(fhash_type__ints_double)          :: vertex_entity_hash, edge_entity_hash, face_entity_hash, domain_entity_hash
type(fhash_type_iterator__ints_double) :: vertex_entity_it, edge_entity_it, face_entity_it, domain_entity_it
type(ints_type)                        :: vertex_entity_key, edge_entity_key, face_entity_key, domain_entity_key
integer                                :: vertex_entity_value, edge_entity_value, face_entity_value, domain_entity_value

type(fhash_type__ints_double)          :: xface1_hash, xface2_hash, yface1_hash, yface2_hash, zface1_hash, zface2_hash
!----------------------------------------------------------------------------------------------------------------------------!
open(unit=12, file=mesh_filename)

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

            do ii = 1, numnp
                read(12,*) (xc(jj,ii), jj = 1, ndm)

                do jj = 1, ndm
                    box_hi(jj) = max(xc(jj,ii), box_hi(jj))
                    box_lo(jj) = min(xc(jj,ii), box_lo(jj))
                enddo
            enddo

            write(iow,*)
            write(6,*)

            write(iow,'(A6,A13,A17,A18)') "dim", "box_length", "box_min", "box_max"
            write(6  ,'(A6,A13,A17,A18)') "dim", "box_length", "box_min", "box_max"
            do jj = 1, ndm
                box_len(jj) = box_hi(jj) - box_lo(jj)
                write(iow,'(I5,2X,3(E16.9,2X))')jj, box_len(jj), box_lo(jj), box_hi(jj)
                write(6  ,'(I5,2X,3(E16.9,2X))')jj, box_len(jj), box_lo(jj), box_hi(jj)
            enddo

            write(6,*)
            write(iow,*)

            box_volume = box_len(1) * box_len(2) * box_len(3)
            write(iow,'(3X,A40,E16.9,A13)')adjl("Box volume:",40),box_volume," [Angstrom^3]"
            write(6  ,'(3X,A40,E16.9,A13)')adjl("Box volume:",40),box_volume," [Angstrom^3]"
        elseif (INDEX(line,"3 vtx # type name")>0) then
            read(12,*)
            read(12,*)
            read(12,*) nen_type_vertex
            read(12,*) numel_type_vertex
            allocate(global_node_id_type_vertex(nen_type_vertex,numel_type_vertex))
            read(12,*)

            do kk = 1, numel_type_vertex
                read(12,*) (global_node_id_type_vertex(pp,kk), pp=1,nen_type_vertex)
            enddo
            read(12,*)
            read(12,'(A100)',IOSTAT=reason) aux_line
            if (INDEX(aux_line," # number of parameter values per element")>0) then
                read(12,*)
                read(12,*)
                read(12,*)
                read(12,*)
            endif
            read(12,*)

            allocate(vertex_entity_key%ints(1))
            call vertex_entity_hash%reserve(numel_type_vertex)
            do kk = 1, numel_type_vertex
                vertex_entity_key%ints(1) = kk
                read(12,*) id_of_entity_it_belongs
                call vertex_entity_hash%get(vertex_entity_key, vertex_entity_value, success)
                if (.NOT.success) call vertex_entity_hash%set(vertex_entity_key, id_of_entity_it_belongs)
            enddo
        elseif ((INDEX(line,"3 edg # type name")>0).OR.(INDEX(line,"4 edg2 # type name")>0)) then
            read(12,*)
            read(12,*)
            read(12,*) nen_type_edge
            read(12,*) numel_type_edge
            allocate(global_node_id_type_edge(nen_type_edge,numel_type_edge))
            read(12,*)

            do kk = 1, numel_type_edge
                read(12,*) (global_node_id_type_edge(pp,kk), pp=1,nen_type_edge)
            enddo
            read(12,*)
            read(12,'(A100)',IOSTAT=reason) aux_line
            if (INDEX(aux_line," # number of parameter values per element")>0) then
                read(12,*)
                read(12,*)
                do ii = 1, numel_type_edge
                    read(12,*)
                enddo
                read(12,*)
                read(12,*)
            endif
            read(12,*)

            allocate(edge_entity_key%ints(1))
            call edge_entity_hash%reserve(numel_type_edge)
            do kk = 1, numel_type_edge
                edge_entity_key%ints(1) = kk
                read(12,*) id_of_entity_it_belongs
                call edge_entity_hash%get(edge_entity_key, edge_entity_value, success)
                if (.NOT.success) call edge_entity_hash%set(edge_entity_key, id_of_entity_it_belongs)
            enddo
        elseif ((INDEX(line,"3 tri # type name")>0).OR.(INDEX(line,"4 tri2 # type name")>0)) then
            read(12,*)
            read(12,*)
            read(12,*) nen_type_face
            read(12,*) numel_type_face
            read(12,*)

            allocate(global_node_id_type_face(nen_type_face,numel_type_face))
            do kk = 1, numel_type_face
                read(12,*) (global_node_id_type_face(pp,kk), pp=1,nen_type_face)
            enddo
            global_node_id_type_face = global_node_id_type_face + 1

            read(12,*)
            read(12,'(A100)',IOSTAT=reason) aux_line
            if ((INDEX(aux_line,"3 # number of parameter values per element")>0).OR.&
               &(INDEX(aux_line,"6 # number of parameter values per element")>0)) then
                read(12,*)
                read(12,*)
                do ii = 1, numel_type_face
                    read(12,*)
                enddo
                read(12,*)
                read(12,*)
            endif
            read(12,*)

            allocate(face_entity_key%ints(1))
            call face_entity_hash%reserve(numel_type_face)
            do kk = 1, numel_type_face
                face_entity_key%ints(1) = kk
                read(12,*) id_of_entity_it_belongs
                call face_entity_hash%get(face_entity_key, face_entity_value, success)
                if (.NOT.success) call face_entity_hash%set(face_entity_key, id_of_entity_it_belongs+1)
            enddo
        elseif ((INDEX(line,"3 tet # type name")>0).OR.(INDEX(line,"4 tet2 # type name")>0)) then
            read(12,*)
            read(12,*)
            read(12,*) nel
            read(12,*) numel
            read(12,*)

            num_of_bulk_pairs = nel * nel * numel

            write(iow,*)
            write(*,*)
            write(iow,'(3X,"Number of mesh points (numnp):         ",I16)') numnp
            write(iow,'(3X,"Number of elements (numel):            ",I16)') numel
            write(iow,'(3X,"Number of nodes per element (nel):     ",I16)') nel
            write(iow,'(3X,"Number of matrix indeces:              ",I16)') num_of_bulk_pairs
            write(6,'(3X,"Number of mesh points (numnp):         ",I16)') numnp
            write(6,'(3X,"Number of elements (numel):            ",I16)') numel
            write(6,'(3X,"Number of nodes per element (nel):     ",I16)') nel
            write(6,'(3X,"Number of matrix indeces:              ",I16)') num_of_bulk_pairs

            allocate(global_node_id_type_domain(nel,numel))
            do ii = 1, numel
                read(12,*) (global_node_id_type_domain(jj,ii), jj = 1, nel)
            enddo
            global_node_id_type_domain = global_node_id_type_domain + 1

            if (nel>4) then
                allocate(temp3(numel))
                do ii = 1, numel
                    temp3(ii) = global_node_id_type_domain(7,ii)
                    global_node_id_type_domain(7,ii) = global_node_id_type_domain(6,ii)
                    global_node_id_type_domain(6,ii) = temp3(ii)
                enddo
            endif

            read(12,*)
            read(12,'(A100)',IOSTAT=reason) aux_line
            if ((INDEX(aux_line,"4 # number of parameter values per element")>0).OR.&
               &(INDEX(aux_line,"10 # number of parameter values per element")>0)) then
                read(12,*)
                read(12,*)
                read(12,*)
                read(12,*)
            endif
            read(12,*)

            allocate(domain_entity_key%ints(1))
            call domain_entity_hash%reserve(numel)
            do kk = 1, numel
                domain_entity_key%ints(1) = kk
                read(12,*) id_of_entity_it_belongs
                call domain_entity_hash%get(domain_entity_key, domain_entity_value, success)
                if (.NOT.success) call domain_entity_hash%set(domain_entity_key, id_of_entity_it_belongs)
            enddo
        endif
    endif
enddo
close(12)

#ifdef DEBUG_OUTPUTS
open(unit=11, file=vertex_elements)
open(unit=22, file=edge_elements)
open(unit=33, file=face_elements)
open(unit=44, file=domain_elements)

call vertex_entity_it%begin(vertex_entity_hash)
call edge_entity_it%begin(edge_entity_hash)
call domain_entity_it%begin(domain_entity_hash)

do kk = 1, vertex_entity_hash%key_count()
    call vertex_entity_it%next(vertex_entity_key, vertex_entity_value)
    write(11,*) vertex_entity_value, vertex_entity_key%ints(1)
enddo

do kk = 1, edge_entity_hash%key_count()
    call edge_entity_it%next(edge_entity_key, edge_entity_value)
    write(22,*) edge_entity_value, edge_entity_key%ints(1)
enddo

do kk = 1, domain_entity_hash%key_count()
    call domain_entity_it%next(domain_entity_key, domain_entity_value)
    write(44,*) domain_entity_value, domain_entity_key%ints(1)
enddo
#endif

if (domain_is_periodic) then
    if (periodic_axis_id(1)) call mesh_face_entities('x', face_entity_hash, xface1_hash, xface2_hash)
    if (periodic_axis_id(2)) call mesh_face_entities('y', face_entity_hash, yface1_hash, yface2_hash)
    if (periodic_axis_id(3)) call mesh_face_entities('z', face_entity_hash, zface1_hash, zface2_hash)

#ifdef DEBUG_OUTPUTS
    call face_entity_it%begin(face_entity_hash)
    do kk = 1, face_entity_hash%key_count()
        call face_entity_it%next(face_entity_key, face_entity_value)
        write(33,*) face_entity_value-1, face_entity_key%ints(1)
    enddo
    close(11)
    close(22)
    close(33)
    close(44)

    if (periodic_axis_id(1)) call mesh_periodic_face_elements('x', xface1_hash, xface2_hash)
    if (periodic_axis_id(2)) call mesh_periodic_face_elements('y', yface1_hash, yface2_hash)
    if (periodic_axis_id(3)) call mesh_periodic_face_elements('z', zface1_hash, zface2_hash)
#endif

    if (periodic_axis_id(1)) call mesh_build_node_pairing(global_node_id_type_face, 'x', xface1_hash, xface2_hash, node_pairing_xx_hash)
    if (periodic_axis_id(2)) call mesh_build_node_pairing(global_node_id_type_face, 'y', yface1_hash, yface2_hash, node_pairing_yy_hash)
    if (periodic_axis_id(3)) call mesh_build_node_pairing(global_node_id_type_face, 'z', zface1_hash, zface2_hash, node_pairing_zz_hash)
endif

allocate(num_of_elems_of_node(numnp))
call mesh_elements_per_node(max_num_of_elems_per_node, num_of_elems_of_node)

num_dest_xx_neighbors = 0
num_dest_yy_neighbors = 0
num_dest_zz_neighbors = 0

if (periodic_axis_id(1)) call mesh_periodic_neighbors(node_pairing_xx_hash, num_dest_xx_neighbors)
if (periodic_axis_id(2)) call mesh_periodic_neighbors(node_pairing_yy_hash, num_dest_yy_neighbors)
if (periodic_axis_id(3)) call mesh_periodic_neighbors(node_pairing_zz_hash, num_dest_zz_neighbors)

total_num_of_node_pairs = num_of_bulk_pairs

if (periodic_axis_id(1)) total_num_of_node_pairs = total_num_of_node_pairs + 2 * (node_pairing_xx_hash%key_count() + num_dest_xx_neighbors)
if (periodic_axis_id(2)) total_num_of_node_pairs = total_num_of_node_pairs + 2 * (node_pairing_yy_hash%key_count() + num_dest_yy_neighbors)
if (periodic_axis_id(3)) total_num_of_node_pairs = total_num_of_node_pairs + 2 * (node_pairing_zz_hash%key_count() + num_dest_zz_neighbors)

allocate(F_m%row(total_num_of_node_pairs))
allocate(F_m%col(total_num_of_node_pairs))
allocate(F_m%g(total_num_of_node_pairs))
allocate(F_m%rh(total_num_of_node_pairs))
allocate(F_m%c(total_num_of_node_pairs))
allocate(F_m%k(total_num_of_node_pairs))
allocate(F_m%w(total_num_of_node_pairs))
allocate(F_m%is_zero(total_num_of_node_pairs))

F_m%g       = 0.d0
F_m%k       = 0.d0
F_m%c       = 0.d0
F_m%rh      = 0.d0
F_m%row     = 0
F_m%col     = 0
F_m%is_zero = .True.

call mesh_bulk_node_pairs(elemcon)

starting_pair = num_of_bulk_pairs + 1
ending_pair   = num_of_bulk_pairs + node_pairing_xx_hash%key_count()
if (periodic_axis_id(1)) call mesh_append_periodic_pairs(elemcon, starting_pair, ending_pair, node_pairing_xx_hash)

do ii = 1, num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count()
     F_m%is_zero(ii) = (node_pair_id(ii)/=ii)
enddo

forward_steps = 0
if (periodic_axis_id(1)) call mesh_append_dest_neighbors(elemcon, forward_steps, num_dest_xx_neighbors, node_pairing_xx_hash)

do ii = num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count() + 1, &
&       num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count() + 2*num_dest_xx_neighbors
     F_m%is_zero(ii) = (node_pair_id(ii)/=ii)
enddo

starting_pair = num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count() + 2*num_dest_xx_neighbors + 1
ending_pair   = num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count() + 2*num_dest_xx_neighbors + node_pairing_yy_hash%key_count()
if (periodic_axis_id(2)) call mesh_append_periodic_pairs(elemcon, starting_pair, ending_pair, node_pairing_yy_hash)

do ii = num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count() + 2*num_dest_xx_neighbors + 1, &
&       num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count() + 2*num_dest_xx_neighbors + 2*node_pairing_yy_hash%key_count()
     F_m%is_zero(ii) = (node_pair_id(ii)/=ii)
enddo

forward_steps = 2*node_pairing_xx_hash%key_count() + 2*num_dest_xx_neighbors
if (periodic_axis_id(2)) call mesh_append_dest_neighbors(elemcon, forward_steps, num_dest_yy_neighbors, node_pairing_yy_hash)

call elemcon%clear()

do ii = num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count() + 2*num_dest_xx_neighbors + 2*node_pairing_yy_hash%key_count() + 1, &
&       num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count() + 2*num_dest_xx_neighbors + 2*node_pairing_yy_hash%key_count() + 2*num_dest_yy_neighbors
     F_m%is_zero(ii) = (node_pair_id(ii)/=ii)
enddo

call mesh_dirichlet_faces(numel_type_face, nen_type_face, global_node_id_type_face, face_entity_hash)

#ifdef DEBUG_OUTPUTS
close(123)
close(13)

open(unit=77, file = com_12)
do ii = 1, total_num_of_node_pairs
     write(77,'(4(2X,I9),2X,L9)') ii, F_m%row(ii), F_m%col(ii), node_pair_id(ii), F_m%is_zero(ii)
enddo
close(77)

open (unit=77, file = inter)
do ii = 1, numel
    write(77,'(11(2X,I9))') (global_node_id_type_domain(jj,ii), jj = 1, nel)
enddo
close(77)

open(77, file = mesh_out)
write(77,'(3(2X,A16))') 'x', 'y', 'z'
do ii = 1, numnp
    write(77,'(3(2X,F16.9))') (xc(jj,ii), jj = 1, ndm)
enddo
close(77)

call mesh_profile()
#endif

#ifdef DEBUG_OUTPUTS
deallocate(vertex_entity_key%ints, edge_entity_key%ints, domain_entity_key%ints)
call vertex_entity_hash%clear()
call edge_entity_hash%clear()
call domain_entity_hash%clear()
#endif

! Deallocate memory
if (nel > 4) deallocate(temp3)
deallocate(global_node_id_type_vertex, global_node_id_type_edge, global_node_id_type_face)

return
!----------------------------------------------------------------------------------------------------------------------------!
end subroutine parser_mesh
