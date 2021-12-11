!--------------------------------------------------------------------------------------------------------------------------------------------!
!                  THIS IS A POST-PROCESSING CODE TO INTERPOLATE THE PROFILES OBTAINED BY THE FEM3D CODE                                     !
!                                 IT APPLIES A SPATIAL THREE-DIMENSIONAL INTERPOLATION USING                                                 !
!                          THE SAME FINITE-ELEMENT SHAPE FUNCTIONS AND MESH WHICH ARE USED IN FEM3D.                                         !
!                                                    DATE: 01/03/2020                                                                        !
!                      AUTHORS: CONSTANTINOS J. REVELAS, ARISTOTELIS P. SGOUROS, APOSTOLIS T. LAKKAS                                         !
!--------------------------------------------------------------------------------------------------------------------------------------------!
program fem3d_spatial_interpolation
!--------------------------------------------------------------------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------------------------------------------------------------------!
integer                              :: i, j, ii, jj, numnp, numel, dummy
integer                              :: nodeId, i_node
integer, parameter                   :: max_el_node=70
integer, allocatable, dimension(:)   :: n_el_node
integer, allocatable, dimension(:,:) :: ix, el_node

real(8)                              :: time_1, time_2, delta_time
real(8)                              :: x_test, y_test, z_test
real(8)                              :: interp, phi_interp
real(8), allocatable, dimension(:)   :: phi_matrix, phi_grafted, wa
real(8), allocatable, dimension(:,:) :: xc
!--------------------------------------------------------------------------------------------------------------------------------------------!
call cpu_time(time_1)

open(unit = 22, file = "in.input.txt")
read(22,*) numnp
read(22,*) numel
close(22)

allocate(xc(3,numnp), ix(4,numel), phi_matrix(numnp), phi_grafted(numnp), wa(numnp))

open (unit = 77, file = "in.elemcon.txt")
read (77,*) ((ix(i,ii), i=1, 4), ii=1, numel)
close(77)

ix = ix + 1

! compute the array elements of node
allocate(n_el_node(numnp))
allocate(el_node(numnp, max_el_node))

n_el_node = 0
do ii = 1, numel
   do jj = 1, 4
      i_node = ix(jj, ii)
      n_el_node(i_node) = n_el_node(i_node) + 1
      el_node(i_node, n_el_node(i_node)) = ii
   enddo
enddo

!do ii = 1, numnp
!   write(123,*) (el_node(ii,jj), jj=1, n_el_node(ii))
!enddo

open(unit = 88, file = "in.profiles.txt")
do i = 1, numnp
    read(88,*) dummy, (xc(j,i), j=1,3), phi_matrix(i), phi_grafted(i), wa(i)
enddo
close(88)

open(unit = 80, file = "o.interp.txt")
!nodeId = 10
!x_test = xc(1,nodeId)
x_test = 0.d0
!y_test = xc(2,nodeId)
y_test = 0.d0
!z_test = xc(3,nodeId)
z_test = 5.d1

do ii = 1, 1000
    z_test     = z_test - 5.d1 / 1000
    phi_interp = interp(numnp, numel, el_node, n_el_node, max_el_node, nodeId, ix, xc, x_test, y_test, z_test, phi_matrix)

    write(80 ,'(I8,4(2X,E16.9))')nodeId, x_test, y_test, z_test, phi_interp
enddo
close(80)

call cpu_time(time_2)

delta_time = time_2 - time_1

write(6, '("cpu time: ", F12.2, "min")') delta_time/6.d01

deallocate(xc, ix, phi_matrix, n_el_node, el_node)
!--------------------------------------------------------------------------------------------------------------------------------------------!
end program fem3d_spatial_interpolation
