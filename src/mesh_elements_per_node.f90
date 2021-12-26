subroutine mesh_elements_per_node(max_num_of_elems_per_node, num_of_elems_of_node)

use parser_vars_mod, only: iow
use geometry_mod, only: numnp, numel, el_node, global_node_id_type_domain

implicit none

integer :: ii, jj, i_node
integer, intent(out) :: max_num_of_elems_per_node
integer, intent(out), dimension(numnp) :: num_of_elems_of_node

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
return
end subroutine mesh_elements_per_node
