real(8) function interp_fem(nodeId, x_interp, y_interp, z_interp, uu)
!--------------------------------------------------------------------------------------------------------------------------------------------!
use geometry
!--------------------------------------------------------------------------------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------------------------------------------------------------------------------!
integer, intent(in)                      :: nodeId
integer                                  :: jj, l, lint, nn, m1, m2, inside, mm, iglob, iloc

real(8), intent(in), dimension(numnp)    :: uu
real(8), intent(in)                      :: x_interp, y_interp, z_interp

real(8), dimension(4,11)                 :: shp
real(8), dimension(5,11)                 :: sv
real(8), dimension(3,4)                  :: xl
real(8), dimension(4)                    :: ul_uu
real(8), dimension(4)                    :: gc, lc
real(8), dimension(4,4)                  :: transf, transf_inv
real(8)                                  :: volel, xsj
real(8)                                  :: u_interp

!--------------------------------------------------------------------------------------------------------------------------------------------!
interface
    pure function matinv4(transf)
    real(8), intent(in), dimension(4,4) :: transf
    real(8),             dimension(4,4) :: matinv4
    end function
end interface
!--------------------------------------------------------------------------------------------------------------------------------------------!

!set up for 11-point quadrature
l = 3
call gausspoints(l, lint, sv)

!loop over elements
do mm = 1, n_el_node(nodeId)
   nn = el_node(nodeId, mm)

    !loop over all nodes of current element
    do iloc = 1, 4
        !find global index of node
        iglob = ix(iloc,nn)

        !copy coordinates from global array to local array concerning current element
        do jj = 1, 3
            xl(jj,iloc) = xc(jj,iglob)
        enddo

        !copy value of solution at current node from global to local array
        ul_uu(iloc)  = uu(iglob)

    enddo

    !initialize accumulator for element volume
    volel = 0.d00

    !loop over all quadrature points in the considered element
    do l = 1, lint
        call tetshp(sv(1,l), xl, ndm, nel, xsj, shp)

        xsj = xsj*sv(5,l)

        volel = volel + xsj
    enddo
    !-------------------------------------------START INTERPOLATION--------------------------------------------------------!
    !check if current plotPoint lies in the specific element
    call isInside(xl, x_interp, y_interp, z_interp, lint, sv, volel, inside, shp)

    if (inside==1) then
        gc(1) = 1.d00
        gc(2) = x_interp
        gc(3) = y_interp
        gc(4) = z_interp

        do m1 = 1, 4
            transf(1,m1) = 1.d00
            transf(2,m1) = xl(1,m1)
            transf(3,m1) = xl(2,m1)
            transf(4,m1) = xl(3,m1)
        enddo

        transf_inv = matinv4(transf)

        do m1 = 1, 4
            lc(m1) = 0.d00
            do m2 = 1, 4
                lc(m1) = lc(m1) + transf_inv(m1,m2) * gc(m2)
            enddo
        enddo

        call tetshp(lc, xl, ndm, nel, xsj, shp)

        u_interp  = 0.d00

        do jj = 1, 4
            u_interp  = u_interp  + shp(4,jj) * ul_uu(jj)
        enddo

        interp_fem = u_interp
        exit

    endif

    !end loop over elements
enddo
!--------------------------------------------------------------------------------------------------------------------------------------------!
end function interp_fem
