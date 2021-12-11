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
use parser_vars_mod,  only: iow, num_of_dirichlet_faces, num_of_nanoparticle_faces,    &
&                       ids_dirichlet_faces, ids_nanopart_faces, periodic_face_id,     &
&                       periodic_axis_id, domain_is_periodic, prof_dim
use geometry_mod,     only: nel, num_of_elems_of_node, box_lo, box_hi, box_len,        &
&                       xc, numnp, numel, el_node, global_node_id_type_domain,         &
&                       is_dirichlet_face, node_belongs_to_dirichlet_face,             &
&                       ndm, node_pair_id, num_of_bulk_pairs, total_num_of_node_pairs, &
&                       node_pairing_xx_hash, node_pairing_yy_hash,                    &
&                       node_pairing_zz_hash, node_pairing_xx_it,                      &
&                       node_pairing_yy_it, node_pairing_zz_it,                        &
&                       node_pairing_xx_key, node_pairing_yy_key,                      &
&                       node_pairing_zz_key, node_pairing_xx_value,                    &
&                       node_pairing_yy_value, node_pairing_zz_value,                  &
&                       num_dest_xx_neighbors
use kcw_mod,          only: F_m
use iofiles_mod,      only: mesh_filename, dir_faces, com_12, inter, mesh_out,         &
&                       mesh_prof, xface1_elements, xface2_elements,                   &
&                       yface1_elements, yface2_elements, zface1_elements,             &
&                       zface2_elements
!----------------------------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------------------------!
character(len=200) :: line, aux_line

integer                            :: node1, node2, elem1, elem2
integer                            :: nbin, ibin, reason, id_of_entity_it_belongs
integer                            :: ii, jj, kk, mm, pp, idummy, i_node
integer, allocatable, dimension(:) :: temp3, prof_1D_node

real(8) :: box_volume = 0.d0, tol = 1.e-8, prof_bin

logical :: success

type(fhash_type__ints_double) :: h
type(ints_type)               :: key
integer                       :: n_keys, key_value

!type(fhash_type__ints_double) :: elemcon
!type(ints_type)             :: elemcon_key
!integer                       :: elemcon_number_of_keys, elemcon_value

type(fhash_type__ints_double)          :: vertex_entity_hash, edge_entity_hash, face_entity_hash, domain_entity_hash
type(fhash_type_iterator__ints_double) :: vertex_entity_it, edge_entity_it, face_entity_it, domain_entity_it
type(ints_type)                        :: vertex_entity_key, edge_entity_key, face_entity_key, domain_entity_key
integer                                :: vertex_entity_value, edge_entity_value, face_entity_value, domain_entity_value

type(fhash_type__ints_double)          :: xface1_hash, xface2_hash, yface1_hash, yface2_hash, zface1_hash, zface2_hash
type(fhash_type_iterator__ints_double) :: xface1_it, xface2_it, yface1_it, yface2_it, zface1_it, zface2_it
type(ints_type)                        :: xface1_key, xface2_key, yface1_key, yface2_key, zface1_key, zface2_key
integer                                :: xface1_value, xface2_value, yface1_value, yface2_value, zface1_value, zface2_value
integer                                :: xface1_size, xface2_size, yface1_size, yface2_size, zface1_size, zface2_size

type(fhash_type__ints_double)          :: node_pairing_xx_hash, node_pairing_yy_hash, node_pairing_zz_hash, node_pairing_all_hash
type(fhash_type_iterator__ints_double) :: node_pairing_xx_it, node_pairing_yy_it, node_pairing_zz_it, node_pairing_all_it
type(ints_type)                        :: node_pairing_xx_key, node_pairing_yy_key, node_pairing_zz_key, node_pairing_all_key
integer                                :: node_pairing_xx_value, node_pairing_yy_value, node_pairing_zz_value, node_pairing_all_value
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

            all_el = nel * nel * numel

            write(iow,*)
            write(*,*)
            write(iow,'(3X,"Number of mesh points (numnp):         ",I16)') numnp
            write(iow,'(3X,"Number of elements (numel):            ",I16)') numel
            write(iow,'(3X,"Number of nodes per element (nel):     ",I16)') nel
            write(iow,'(3X,"Number of matrix indeces:              ",I16)') all_el
            write(6,'(3X,"Number of mesh points (numnp):         ",I16)') numnp
            write(6,'(3X,"Number of elements (numel):            ",I16)') numel
            write(6,'(3X,"Number of nodes per element (nel):     ",I16)') nel
            write(6,'(3X,"Number of matrix indeces:              ",I16)') all_el

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
    ! Calculate the size of the face entity hashes
    xface1_size = 0
    xface2_size = 0
    yface1_size = 0
    yface2_size = 0
    zface1_size = 0
    zface2_size = 0

    call face_entity_it%begin(face_entity_hash)
    do kk = 1, face_entity_hash%key_count()
        call face_entity_it%next(face_entity_key, face_entity_value)
#ifdef DEBUG_OUTPUTS
        write(33,*) face_entity_value-1, face_entity_key%ints(1)
#endif
        if (periodic_axis_id(1)) then
            if (face_entity_value==periodic_face_id(1)) xface1_size = xface1_size + 1
            if (face_entity_value==periodic_face_id(2)) xface2_size = xface2_size + 1
        endif
        if (periodic_axis_id(2)) then
            if (face_entity_value==periodic_face_id(3)) yface1_size = yface1_size + 1
            if (face_entity_value==periodic_face_id(4)) yface2_size = yface2_size + 1
        endif
        if (periodic_axis_id(3)) then
            if (face_entity_value==periodic_face_id(5)) zface1_size = zface1_size + 1
            if (face_entity_value==periodic_face_id(6)) zface2_size = zface2_size + 1
        endif
    enddo

#ifdef DEBUG_OUTPUTS
    close(11)
    close(22)
    close(33)
    close(44)
#endif

    if (periodic_axis_id(1)) then
        allocate(xface1_key%ints(1))
        allocate(xface2_key%ints(1))
        call xface1_hash%reserve(xface1_size)
        call xface2_hash%reserve(xface2_size)
    endif
    if (periodic_axis_id(2)) then
        allocate(yface1_key%ints(1))
        allocate(yface2_key%ints(1))
        call yface1_hash%reserve(yface1_size)
        call yface2_hash%reserve(yface2_size)
    endif
    if (periodic_axis_id(3)) then
        allocate(zface1_key%ints(1))
        allocate(zface2_key%ints(1))
        call zface1_hash%reserve(zface1_size)
        call zface2_hash%reserve(zface2_size)
    endif

    ! Fill the face entity hashes
    call face_entity_it%begin(face_entity_hash)
    do kk = 1, face_entity_hash%key_count()
        call face_entity_it%next(face_entity_key, face_entity_value)
        if (face_entity_value==periodic_face_id(1)) then
            xface1_key%ints(1) = face_entity_key%ints(1)
            call xface1_hash%set(xface1_key, face_entity_value)
            cycle
        endif
        if (face_entity_value==periodic_face_id(2)) then
            xface2_key%ints(1) = face_entity_key%ints(1)
            call xface2_hash%set(xface2_key, face_entity_value)
            cycle
        endif
        if (face_entity_value==periodic_face_id(3)) then
            yface1_key%ints(1) = face_entity_key%ints(1)
            call yface1_hash%set(yface1_key, face_entity_value)
            cycle
        endif
        if (face_entity_value==periodic_face_id(4)) then
            yface2_key%ints(1) = face_entity_key%ints(1)
            call yface2_hash%set(yface2_key, face_entity_value)
            cycle
        endif
        if (face_entity_value==periodic_face_id(5)) then
            zface1_key%ints(1) = face_entity_key%ints(1)
            call zface1_hash%set(zface1_key, face_entity_value)
            cycle
        endif
        if (face_entity_value==periodic_face_id(6)) then
            zface2_key%ints(1) = face_entity_key%ints(1)
            call zface2_hash%set(zface2_key, face_entity_value)
            cycle
        endif
    enddo

#ifdef DEBUG_OUTPUTS
    if (periodic_axis_id(1)) then
        open(unit=111, file=xface1_elements)
        open(unit=222, file=xface2_elements)

        call xface1_it%begin(xface1_hash)
        call xface2_it%begin(xface2_hash)

        do kk = 1, xface1_hash%key_count()
            call xface1_it%next(xface1_key, xface1_value)
            write(111,*) xface1_value, xface1_key%ints(1)
        enddo
        do kk = 1, xface2_hash%key_count()
            call xface2_it%next(xface2_key, xface2_value)
            write(222,*) xface2_value, xface2_key%ints(1)
        enddo

        close(111)
        close(222)
    endif

    if (periodic_axis_id(2)) then
        open(unit=333, file=yface1_elements)
        open(unit=444, file=yface2_elements)

        call yface1_it%begin(yface1_hash)
        call yface2_it%begin(yface2_hash)

        do kk = 1, yface1_hash%key_count()
            call yface1_it%next(yface1_key, yface1_value)
            write(333,*) yface1_value, yface1_key%ints(1)
        enddo
        do kk = 1, yface2_hash%key_count()
            call yface2_it%next(yface2_key, yface2_value)
            write(444,*) yface2_value, yface2_key%ints(1)
        enddo

        close(333)
        close(444)
    endif

    if (periodic_axis_id(3)) then
        open(unit=555, file=zface1_elements)
        open(unit=666, file=zface2_elements)

        call zface1_it%begin(zface1_hash)
        call zface2_it%begin(zface2_hash)

        do kk = 1, zface1_hash%key_count()
            call zface1_it%next(zface1_key, zface1_value)
            write(555,*) zface1_value, zface1_key%ints(1)
        enddo
        do kk = 1, zface2_hash%key_count()
            call zface2_it%next(zface2_key, zface2_value)
            write(666,*) zface2_value, zface2_key%ints(1)
        enddo

        close(555)
        close(666)
    endif
#endif

    allocate(node_pairing_all_key%ints(1))
    call node_pairing_all_hash%reserve(xface1_size+yface1_size+zface1_size)

    ! Build xx pairing
    if (periodic_axis_id(1)) then
        allocate(node_pairing_xx_key%ints(1))
        call node_pairing_xx_hash%reserve(xface1_hash%key_count())
        call xface1_it%begin(xface1_hash)
        do kk = 1, xface1_hash%key_count()
            call xface1_it%next(xface1_key, xface1_value)
            elem1 = xface1_key%ints(1)

            call xface2_it%begin(xface2_hash)
            do mm = 1, xface2_hash%key_count()
                call xface2_it%next(xface2_key, xface2_value)
                elem2 = xface2_key%ints(1)
                do ii = 1, 3
                    node1 = global_node_id_type_face(ii,elem1)
                    do jj = 1, 3
                        node2 = global_node_id_type_face(jj, elem2)
                        if ((ABS(xc(2,node1)-xc(2,node2))<1e-12).AND.(ABS(xc(3,node1)-xc(3,node2))<1e-12)) then
                            node_pairing_xx_key%ints(1) = node1
                            call node_pairing_xx_hash%set(node_pairing_xx_key, node2)
                            node_pairing_xx_key%ints(1) = node2
                            call node_pairing_xx_hash%set(node_pairing_xx_key, node1)

                            node_pairing_all_key%ints(1) = node1
                            call node_pairing_all_hash%set(node_pairing_all_key, node2)
                            node_pairing_all_key%ints(1) = node2
                            call node_pairing_all_hash%set(node_pairing_all_key, node1)
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif

    ! Build yy pairing
    if (periodic_axis_id(2)) then
        allocate(node_pairing_yy_key%ints(1))
        call node_pairing_yy_hash%reserve(yface1_hash%key_count())
        call yface1_it%begin(yface1_hash)
        do kk = 1, yface1_hash%key_count()
            call yface1_it%next(yface1_key, yface1_value)
            elem1 = yface1_key%ints(1)

            call yface2_it%begin(yface2_hash)
            do mm = 1, yface2_hash%key_count()
                call yface2_it%next(yface2_key, yface2_value)
                elem2 = yface2_key%ints(1)
                do ii = 1, 3
                    node1 = global_node_id_type_face(ii,elem1)
                    do jj = 1, 3
                        node2 = global_node_id_type_face(jj, elem2)
                        if ((ABS(xc(1,node1)-xc(1,node2))<1e-12).AND.(ABS(xc(3,node1)-xc(3,node2))<1e-12)) then
                            node_pairing_yy_key%ints(1) = node1
                            call node_pairing_yy_hash%set(node_pairing_yy_key, node2)
                            node_pairing_yy_key%ints(1) = node2
                            call node_pairing_yy_hash%set(node_pairing_yy_key, node1)

                            node_pairing_all_key%ints(1) = node1
                            call node_pairing_all_hash%set(node_pairing_all_key, node2)
                            node_pairing_all_key%ints(1) = node2
                            call node_pairing_all_hash%set(node_pairing_all_key, node1)
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif

    ! Build zz pairing
    if (periodic_axis_id(3)) then
        allocate(node_pairing_zz_key%ints(1))
        call node_pairing_zz_hash%reserve(zface1_hash%key_count())
        call zface1_it%begin(zface1_hash)
        do kk = 1, zface1_hash%key_count()
            call zface1_it%next(zface1_key, zface1_value)
            elem1 = zface1_key%ints(1)

            call zface2_it%begin(zface2_hash)
            do mm = 1, zface2_hash%key_count()
                call zface2_it%next(zface2_key, zface2_value)
                elem2 = zface2_key%ints(1)
                do ii = 1, 3
                    node1 = global_node_id_type_face(ii,elem1)
                    do jj = 1, 3
                        node2 = global_node_id_type_face(jj, elem2)
                        if ((ABS(xc(1,node1)-xc(1,node2))<1e-12).AND.(ABS(xc(2,node1)-xc(2,node2))<1e-12)) then
                            node_pairing_zz_key%ints(1) = node1
                            call node_pairing_zz_hash%set(node_pairing_zz_key, node2)
                            node_pairing_zz_key%ints(1) = node2
                            call node_pairing_zz_hash%set(node_pairing_zz_key, node1)

                            node_pairing_all_key%ints(1) = node1
                            call node_pairing_all_hash%set(node_pairing_all_key, node2)
                            node_pairing_all_key%ints(1) = node2
                            call node_pairing_all_hash%set(node_pairing_all_key, node1)
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif

#ifdef DEBUG_OUTPUTS
    if (periodic_axis_id(1)) then
        open(unit=1111,file="d.node_pairing_xx")

        call node_pairing_xx_it%begin(node_pairing_xx_hash)

        do kk = 1, node_pairing_xx_hash%key_count()
            call node_pairing_xx_it%next(node_pairing_xx_key, node_pairing_xx_value)
            write(1111,*) node_pairing_xx_key%ints(1), node_pairing_xx_value
        enddo

        close(1111)
    endif
    if (periodic_axis_id(2)) then
        open(unit=2222,file="d.node_pairing_yy")

        call node_pairing_yy_it%begin(node_pairing_yy_hash)

        do kk = 1, node_pairing_yy_hash%key_count()
            call node_pairing_yy_it%next(node_pairing_yy_key, node_pairing_yy_value)
            write(2222,*) node_pairing_yy_key%ints(1), node_pairing_yy_value
        enddo

        close(2222)
    endif
    if (periodic_axis_id(3)) then
        open(unit=3333,file="d.node_pairing_zz")

        call node_pairing_zz_it%begin(node_pairing_zz_hash)

        do kk = 1, node_pairing_zz_hash%key_count()
            call node_pairing_zz_it%next(node_pairing_zz_key, node_pairing_zz_value)
            write(3333,*) node_pairing_zz_key%ints(1), node_pairing_zz_value
        enddo

        close(3333)
    endif

    call node_pairing_all_it%begin(node_pairing_all_hash)
    do kk = 1, node_pairing_all_hash%key_count()
        call node_pairing_all_it%next(node_pairing_all_key, node_pairing_all_value)
        write(999,*) node_pairing_all_key%ints(1), node_pairing_all_value
    enddo
#endif
endif !domain_is_periodic

! Compute the maximum number of elements per node
allocate(n_el_node(numnp))
n_el_node = 0

max_el_node = 0
do ii = 1, numel
   do jj = 1, 4
      i_node = global_node_id_type_domain(jj, ii)
      n_el_node(i_node) = n_el_node(i_node) + 1
      max_el_node = MAX(max_el_node, n_el_node(i_node))
   enddo
enddo

write(iow,'(3X,"Maximum number of elements per node: ",I18)') max_el_node
write(6  ,'(3X,"Maximum number of elements per node: ",I18)') max_el_node

allocate(el_node(numnp, max_el_node))
el_node = 0

n_el_node = 0
do ii = 1, numel
   do jj = 1, 4
      i_node = global_node_id_type_domain(jj, ii)
      n_el_node(i_node) = n_el_node(i_node) + 1
      el_node(i_node, n_el_node(i_node)) = ii
   enddo
enddo

! Allocate and initialize arrays for matrix assembly
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

! Assembly the con_12 (hash) matrix
allocate(key%ints(2)) ! Each key is defined by a pair (2) of nodes

! Total number of required keys
n_keys = nel * numel
call h%reserve(n_keys*2)

all_el = 0
do mm = 1, numel
    do jj = 1, nel
        do ii = 1, nel
            all_el = all_el + 1

            F_m%row(all_el) = global_node_id_type_domain(jj,mm)
            F_m%col(all_el) = global_node_id_type_domain(ii,mm)

            ! Define the pair of nodes to be examined and assigned a key_value
            key%ints(1) = global_node_id_type_domain(jj,mm)
            key%ints(2) = global_node_id_type_domain(ii,mm)

            ! Assign key_value to the pair
            call h%get(key, key_value, success)

            if (success) then
               con_l2(all_el) = key_value  ! This pair has already been met, thus assigned a key_value
            else
               call h%set(key, all_el)     ! Store the new key_value for next iteration's check
               con_l2(all_el) = all_el     ! This pair is met for the first time
            endif
        enddo
    enddo
enddo
call h%clear()

! Determine all elements belonging to dirichlet faces
allocate(node_in_q0_face(numnp))
node_in_q0_face = .false.

#ifdef DEBUG_OUTPUTS
open(unit=123, file = dir_faces)
#endif

is_dir_face = .false.

do jj = 1, numel_type_face
    face_entity_key%ints(1) = jj
    call face_entity_hash%get(face_entity_key, face_entity_value)
    do ii = 1, n_dirichlet_faces
        if (face_entity_value==ids_dirichlet_faces(ii)) then
            do pp = 1, nen_type_face
                idummy = global_node_id_type_face(pp,jj)

                ! If the node belongs to a dirichlet face
                node_in_q0_face(idummy) = .True.

                ! Find if a node is located at a corner
                kk = 0
                do mm = 1, ndm
                    if (DABS(xc(mm, idummy) - box_lo(mm)) < tol) then
                        kk = kk + 1
                    endif
                    if (DABS(xc(mm, idummy) - box_hi(mm)) < tol) then
                        kk = kk + 1
                    endif
                enddo

                ! If a node is located at a corner skip the loop
                if (kk > 1) cycle

                do mm = 1, ndm
                    if (DABS(xc(mm, idummy) - box_lo(mm)) < tol) then
                        is_dir_face(mm,1) = .true.
                    endif
                    if (DABS(xc(mm, idummy) - box_hi(mm)) < tol) then
                        is_dir_face(mm,2) = .true.
                    endif
                enddo
#ifdef DEBUG_OUTPUTS
                write(123,'(3(E16.8),6(L3))') xc(1, idummy), xc(2, idummy), xc(3, idummy), &
    &                                 (is_dir_face(kk,1), kk = 1,3), (is_dir_face(kk,2), kk = 1,3)
#endif
            enddo
        endif
    enddo

    ! Nanoparticles section
    do ii = 1, n_nanopart_faces
        if (face_entity_value==ids_nanopart_faces(ii)) then
            do kk = 1, nen_type_face
                idummy = global_node_id_type_face(kk,jj)

                ! If the node belongs to a dirichlet face
                node_in_q0_face(idummy) = .True.
            enddo
        endif
    enddo
enddo

#ifdef DEBUG_OUTPUTS
close(123)
close(13)

open(unit=77, file = com_12)
do ii = 1, all_el
     write(77,'(4(2X,I9),2X,L9)') ii, F_m%row(ii), F_m%col(ii), con_l2(ii), con_l2(ii)/=ii
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

!TODO: This should be encapsulated into a subroutine
prof_bin = 0.5d0
nbin = NINT((box_hi(prof_dim) - box_lo(prof_dim)) / prof_bin) + 1

allocate(prof_1D_node(nbin))

prof_1D_node=0
do ii = 1, numnp
    ibin = NINT((xc(prof_dim,ii) - box_lo(prof_dim))/prof_bin) + 1
    prof_1D_node(ibin) = prof_1D_node(ibin) + 1
enddo

open(77, file = mesh_prof)
do ii = 1, nbin
    write(77,*) ii, prof_1D_node(ii)
enddo
close(77)
#endif

#ifdef DEBUG_OUTPUTS
deallocate(vertex_entity_key%ints)
deallocate(edge_entity_key%ints)
deallocate(domain_entity_key%ints)

call vertex_entity_hash%clear()
call edge_entity_hash%clear()
call domain_entity_hash%clear()

if (domain_is_periodic) then
    if (periodic_axis_id(1)) then
        deallocate(xface1_key%ints)
        deallocate(xface2_key%ints)
        call xface1_hash%clear()
        call xface2_hash%clear()
    endif
    if (periodic_axis_id(2)) then
        deallocate(yface1_key%ints)
        deallocate(yface2_key%ints)
        call yface1_hash%clear()
        call yface2_hash%clear()
    endif
    if (periodic_axis_id(3)) then
        deallocate(zface1_key%ints)
        deallocate(zface2_key%ints)
        call zface1_hash%clear()
        call zface2_hash%clear()
    endif
endif
#endif

return
!----------------------------------------------------------------------------------------------------------------------------!
end subroutine parser_mesh
