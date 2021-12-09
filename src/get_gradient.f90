real(8) function get_gradient(node_id, phi)
!----------------------------------------------------------------------------------------------------------!
use geometry
use constants, only: A_to_m
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: node_id 

real(8), intent(in), dimension(numnp) :: phi
real(8)                               :: fem_interpolation
real(8)                               :: phi_plus_dx = 0.d0, phi_minus_dx = 0.d0
real(8)                               :: phi_plus_dy = 0.d0, phi_minus_dy = 0.d0
real(8)                               :: phi_plus_dz = 0.d0, phi_minus_dz = 0.d0   
real(8)                               :: x_node_id = 0.d0, y_node_id = 0.d0, z_node_id = 0.d0
real(8)                               :: dphi2_dx2 = 0.d0, dphi2_dy2 = 0.d0, dphi2_dz2 = 0.d0
!----------------------------------------------------------------------------------------------------------!
x_node_id = xc(1,node_id)
y_node_id = xc(2,node_id)
z_node_id = xc(3,node_id)

phi_plus_dx  = fem_interpolation(node_id, x_node_id+dx, y_node_id, z_node_id, phi)
phi_minus_dx = fem_interpolation(node_id, x_node_id-dx, y_node_id, z_node_id, phi)

phi_plus_dy  = fem_interpolation(node_id, x_node_id, y_node_id+dy, z_node_id, phi)
phi_minus_dy = fem_interpolation(node_id, x_node_id, y_node_id-dy, z_node_id, phi)

phi_plus_dz  = fem_interpolation(node_id, x_node_id, y_node_id, z_node_id+dz, phi)
phi_minus_dz = fem_interpolation(node_id, x_node_id, y_node_id, z_node_id-dz, phi)

dphi2_dx2 = (phi_plus_dx - 2.d00 * phi(node_id) + phi_minus_dx) / (dx*A_to_m)**2.d0
dphi2_dy2 = (phi_plus_dy - 2.d00 * phi(node_id) + phi_minus_dy) / (dy*A_to_m)**2.d0
dphi2_dz2 = (phi_plus_dz - 2.d00 * phi(node_id) + phi_minus_dz) / (dz*A_to_m)**2.d0

get_gradient = dphi2_dx2 + dphi2_dy2 + dphi2_dz2

!write(123,'(2(E16.5,2X))') z_node_id, phi(node_id)!, phi_plus_dx, phi_minus_dx, phi_plus_dy, phi_minus_dy, phi_plus_dz, phi_minus_dz, &
           !& dphi2_dx2, dphi2_dy2, dphi2_dz2, compute_gradient
return
!----------------------------------------------------------------------------------------------------------!
end function get_gradient
