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
use parser_vars_mod,  only: iow, periodic_face_id, periodic_axis_id, domain_is_periodic
use geometry_mod,     only: nel, num_of_elems_of_node, box_lo, box_hi, box_len,        &
&                       xc, numnp, numel, el_node, global_node_id_type_domain,         &
&                       ndm, node_pair_id, num_of_bulk_pairs, total_num_of_node_pairs, &
&                       node_pairing_xx_hash, node_pairing_yy_hash,                    &
&                       node_pairing_zz_hash, node_pairing_xx_it,                      &
&                       node_pairing_yy_it, node_pairing_zz_it,                        &
&                       node_pairing_xx_key, node_pairing_yy_key,                      &
&                       node_pairing_zz_key, node_pairing_xx_value,                    &
&                       node_pairing_yy_value, node_pairing_zz_value,                  &
&                       num_dest_xx_neighbors, num_dest_yy_neighbors,                  &
&                       num_dest_zz_neighbors
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

integer                              :: node1, node2, elem1, elem2, i_node
integer                              :: ii, jj, kk, mm, pp
integer                              :: aux1, aux2, aux3, aux4
integer                              :: reason, id_of_entity_it_belongs, max_num_of_elems_per_node
integer                              :: nen_type_vertex, numel_type_vertex
integer                              :: nen_type_edge, numel_type_edge
integer                              :: nen_type_face, numel_type_face
integer                              :: source_xx, dest_xx, node_pair
integer                              :: source_yy, dest_yy
integer, allocatable, dimension(:)   :: temp3
integer, allocatable, dimension(:,:) :: global_node_id_type_vertex, global_node_id_type_edge, global_node_id_type_face

real(8) :: box_volume = 0.d0

logical :: success

type(fhash_type__ints_double) :: elemcon
type(ints_type)               :: elemcon_key
integer                       :: elemcon_value

type(fhash_type__ints_double)          :: vertex_entity_hash, edge_entity_hash, face_entity_hash, domain_entity_hash
type(fhash_type_iterator__ints_double) :: vertex_entity_it, edge_entity_it, face_entity_it, domain_entity_it
type(ints_type)                        :: vertex_entity_key, edge_entity_key, face_entity_key, domain_entity_key
integer                                :: vertex_entity_value, edge_entity_value, face_entity_value, domain_entity_value

type(fhash_type__ints_double)          :: xface1_hash, xface2_hash, yface1_hash, yface2_hash, zface1_hash, zface2_hash
type(fhash_type_iterator__ints_double) :: xface1_it, xface2_it, yface1_it, yface2_it, zface1_it, zface2_it
type(ints_type)                        :: xface1_key, xface2_key, yface1_key, yface2_key, zface1_key, zface2_key
integer                                :: xface1_value, xface2_value, yface1_value, yface2_value, zface1_value, zface2_value
integer                                :: xface1_size, xface2_size, yface1_size, yface2_size, zface1_size, zface2_size
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
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif

#ifdef DEBUG_OUTPUTS
    if (periodic_axis_id(1)) then
        open(unit=1111, file=node_pairing_xx)

        call node_pairing_xx_it%begin(node_pairing_xx_hash)

        do kk = 1, node_pairing_xx_hash%key_count()
            call node_pairing_xx_it%next(node_pairing_xx_key, node_pairing_xx_value)
            write(1111,*) node_pairing_xx_key%ints(1), node_pairing_xx_value
        enddo

        close(1111)
    endif
    if (periodic_axis_id(2)) then
        open(unit=2222, file=node_pairing_yy)

        call node_pairing_yy_it%begin(node_pairing_yy_hash)

        do kk = 1, node_pairing_yy_hash%key_count()
            call node_pairing_yy_it%next(node_pairing_yy_key, node_pairing_yy_value)
            write(2222,*) node_pairing_yy_key%ints(1), node_pairing_yy_value
        enddo

        close(2222)
    endif
    if (periodic_axis_id(3)) then
        open(unit=3333, file=node_pairing_zz)

        call node_pairing_zz_it%begin(node_pairing_zz_hash)

        do kk = 1, node_pairing_zz_hash%key_count()
            call node_pairing_zz_it%next(node_pairing_zz_key, node_pairing_zz_value)
            write(3333,*) node_pairing_zz_key%ints(1), node_pairing_zz_value
        enddo

        close(3333)
    endif
#endif
endif !domain_is_periodic

! Compute the maximum number of elements per node
allocate(num_of_elems_of_node(numnp))
num_of_elems_of_node = 0

max_num_of_elems_per_node = 0
do ii = 1, numel
   do jj = 1, 4
      i_node = global_node_id_type_domain(jj, ii)

      num_of_elems_of_node(i_node) = num_of_elems_of_node(i_node) + 1
      max_num_of_elems_per_node    = MAX(max_num_of_elems_per_node, num_of_elems_of_node(i_node))
   enddo
enddo

write(iow,'(3X,"Maximum number of elements per node: ",I18)') max_num_of_elems_per_node
write(6  ,'(3X,"Maximum number of elements per node: ",I18)') max_num_of_elems_per_node

allocate(el_node(numnp, max_num_of_elems_per_node))
el_node = 0

num_of_elems_of_node = 0
do ii = 1, numel
   do jj = 1, 4
      i_node = global_node_id_type_domain(jj, ii)

      num_of_elems_of_node(i_node)                  = num_of_elems_of_node(i_node) + 1
      el_node(i_node, num_of_elems_of_node(i_node)) = ii
   enddo
enddo

num_dest_xx_neighbors = 0
if (periodic_axis_id(1)) then
    call node_pairing_xx_it%begin(node_pairing_xx_hash)

    do kk = 1, node_pairing_xx_hash%key_count()
        call node_pairing_xx_it%next(node_pairing_xx_key, node_pairing_xx_value)

        source_xx = node_pairing_xx_key%ints(1)
        dest_xx   = node_pairing_xx_value

        num_dest_xx_neighbors = num_dest_xx_neighbors + 3 * num_of_elems_of_node(dest_xx)
    enddo
endif

num_dest_yy_neighbors = 0
if (periodic_axis_id(2)) then
    call node_pairing_yy_it%begin(node_pairing_yy_hash)

    do kk = 1, node_pairing_yy_hash%key_count()
        call node_pairing_yy_it%next(node_pairing_yy_key, node_pairing_yy_value)

        source_yy = node_pairing_yy_key%ints(1)
        dest_yy   = node_pairing_yy_value

        num_dest_yy_neighbors = num_dest_yy_neighbors + 3 * num_of_elems_of_node(dest_yy)
    enddo
endif

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

allocate(elemcon_key%ints(2))

! Append xx periodic pairs in F_m struct
if (periodic_axis_id(1)) then
    call node_pairing_xx_it%begin(node_pairing_xx_hash)

    do kk = num_of_bulk_pairs + 1, num_of_bulk_pairs + node_pairing_xx_hash%key_count()
        call node_pairing_xx_it%next(node_pairing_xx_key, node_pairing_xx_value)

        source_xx = node_pairing_xx_key%ints(1)
        dest_xx   = node_pairing_xx_value

        ! Append pair
        F_m%row(kk) = source_xx
        F_m%col(kk) = dest_xx

        elemcon_key%ints(1) = source_xx
        elemcon_key%ints(2) = dest_xx

        call elemcon%get(elemcon_key, elemcon_value, success)

        if (success) then
            node_pair_id(kk) = elemcon_value
        else
            call elemcon%set(elemcon_key, kk)
            node_pair_id(kk) = kk
        endif

        ! Append inverse pair
        pp = kk + node_pairing_xx_hash%key_count()

        F_m%row(pp) = dest_xx
        F_m%col(pp) = source_xx

        elemcon_key%ints(1) = dest_xx
        elemcon_key%ints(2) = source_xx

        call elemcon%get(elemcon_key, elemcon_value, success)

        if (success) then
            node_pair_id(pp) = elemcon_value
        else
            call elemcon%set(elemcon_key, pp)
            node_pair_id(pp) = pp
        endif
    enddo
endif

aux1 = num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count()

do ii = 1, aux1
     F_m%is_zero(ii) = (node_pair_id(ii)/=ii)
enddo

! Append dest_xx neighbors in F_m struct
if (periodic_axis_id(1)) then
    call node_pairing_xx_it%begin(node_pairing_xx_hash)

    do kk = 1, node_pairing_xx_hash%key_count()
        call node_pairing_xx_it%next(node_pairing_xx_key, node_pairing_xx_value)

        source_xx = node_pairing_xx_key%ints(1)
        dest_xx   = node_pairing_xx_value

        do node_pair = 1, num_of_bulk_pairs
            if ((F_m%col(node_pair)==dest_xx).and.(F_m%row(node_pair).ne.dest_xx)) then

                if (F_m%is_zero(node_pair)) cycle

                aux2 = aux2 + 1

                pp = aux1 + aux2

                ! Append pair
                F_m%row(pp) = source_xx
                F_m%col(pp) = F_m%row(node_pair)

                elemcon_key%ints(1) = source_xx
                elemcon_key%ints(2) = F_m%row(node_pair)

                call elemcon%get(elemcon_key, elemcon_value, success)

                if (success) then
                    node_pair_id(pp) = elemcon_value
                else
                    call elemcon%set(elemcon_key, pp)
                    node_pair_id(pp) = pp
                endif

                ! Append inverse pair
                F_m%row(pp + num_dest_xx_neighbors) = F_m%row(node_pair)
                F_m%col(pp + num_dest_xx_neighbors) = source_xx

                elemcon_key%ints(1) = F_m%row(node_pair)
                elemcon_key%ints(2) = source_xx

                call elemcon%get(elemcon_key, elemcon_value, success)

                if (success) then
                    node_pair_id(pp + num_dest_xx_neighbors) = elemcon_value
                else
                    call elemcon%set(elemcon_key, pp + num_dest_xx_neighbors)
                    node_pair_id(pp + num_dest_xx_neighbors) = pp + num_dest_xx_neighbors
                endif
            endif
        enddo
    enddo
endif

do ii = aux1 + 1, aux1 + 2*num_dest_xx_neighbors
     F_m%is_zero(ii) = (node_pair_id(ii)/=ii)
enddo

aux3 = aux1 + 2*num_dest_xx_neighbors

! Append yy periodic pairs in F_m struct
if (periodic_axis_id(2)) then
    call node_pairing_yy_it%begin(node_pairing_yy_hash)

    do kk = aux3 + 1, aux3 + node_pairing_yy_hash%key_count()
        call node_pairing_yy_it%next(node_pairing_yy_key, node_pairing_yy_value)

        source_yy = node_pairing_yy_key%ints(1)
        dest_yy   = node_pairing_yy_value

        ! Append pair
        F_m%row(kk) = source_yy
        F_m%col(kk) = dest_yy

        elemcon_key%ints(1) = source_yy
        elemcon_key%ints(2) = dest_yy

        call elemcon%get(elemcon_key, elemcon_value, success)

        if (success) then
            node_pair_id(kk) = elemcon_value
        else
            call elemcon%set(elemcon_key, kk)
            node_pair_id(kk) = kk
        endif

        ! Append inverse pair
        pp = kk + node_pairing_yy_hash%key_count()

        F_m%row(pp) = dest_yy
        F_m%col(pp) = source_yy

        elemcon_key%ints(1) = dest_yy
        elemcon_key%ints(2) = source_yy

        call elemcon%get(elemcon_key, elemcon_value, success)

        if (success) then
            node_pair_id(pp) = elemcon_value
        else
            call elemcon%set(elemcon_key, pp)
            node_pair_id(pp) = pp
        endif
    enddo
endif

aux4 = aux3 + 2*node_pairing_yy_hash%key_count()

do ii = aux3+1, aux4
     F_m%is_zero(ii) = (node_pair_id(ii)/=ii)
enddo

aux2 = 0

!Append dest_yy neighbors in F_m struct
if (periodic_axis_id(2)) then
    call node_pairing_yy_it%begin(node_pairing_yy_hash)

    do kk = 1, node_pairing_yy_hash%key_count()
        call node_pairing_yy_it%next(node_pairing_yy_key, node_pairing_yy_value)

        source_yy = node_pairing_yy_key%ints(1)
        dest_yy   = node_pairing_yy_value

        do node_pair = 1, num_of_bulk_pairs + 2*node_pairing_xx_hash%key_count() + 2*num_dest_xx_neighbors
            if ((F_m%col(node_pair)==dest_yy).and.(F_m%row(node_pair).ne.dest_yy)) then

                if (F_m%is_zero(node_pair)) cycle

                aux2 = aux2 + 1

                pp = aux4 + aux2

                ! Append pair
                F_m%row(pp) = source_yy
                F_m%col(pp) = F_m%row(node_pair)

                elemcon_key%ints(1) = source_yy
                elemcon_key%ints(2) = F_m%row(node_pair)

                call elemcon%get(elemcon_key, elemcon_value, success)

                if (success) then
                    node_pair_id(pp) = elemcon_value
                else
                    call elemcon%set(elemcon_key, pp)
                    node_pair_id(pp) = pp
                endif

                ! Append inverse pair
                F_m%row(pp + num_dest_yy_neighbors) = F_m%row(node_pair)
                F_m%col(pp + num_dest_yy_neighbors) = source_yy

                elemcon_key%ints(1) = F_m%row(node_pair)
                elemcon_key%ints(2) = source_yy

                call elemcon%get(elemcon_key, elemcon_value, success)

                if (success) then
                    node_pair_id(pp + num_dest_yy_neighbors) = elemcon_value
                else
                    call elemcon%set(elemcon_key, pp + num_dest_yy_neighbors)
                    node_pair_id(pp + num_dest_yy_neighbors) = pp + num_dest_yy_neighbors
                endif
            endif
        enddo
    enddo
endif
call elemcon%clear()

do ii = aux4 + 1, aux4 + 2*num_dest_yy_neighbors
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

deallocate(elemcon_key%ints)
return
!----------------------------------------------------------------------------------------------------------------------------!
end subroutine parser_mesh
