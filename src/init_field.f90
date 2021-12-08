subroutine init_field(Ufield, wa)
!------------------------------------------------------------------------------------------------------!
use parser_vars
use geometry
use error_handing
use write_helper
use force_fields
use iofiles
use constants
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: k1, m1, n1

real(8), intent(out), dimension(numnp) :: Ufield, wa
real(8)                                :: number_density = 0.d0
real(8)                                :: r_centers = 0.d0, r_center_surf = 0.d0, r_surf = 0.d0
real(8)                                :: radius_np_actual = 0.d0, radius_pol = 0.d0
real(8)                                :: Urep = 0.d0, Uatt = 0.d0
!------------------------------------------------------------------------------------------------------!
open(unit=211, file = usolid)
write(211,'(3A16,2A19)') "r_center_surf", "r_surf", "Uatt", "Urep", "Utot"

Ufield = 0.d0

do k1 = 1, numnp
    !loop over all dirichlet faces
    do m1 = 1, 3
        do n1 = 1, 2
            if (is_dir_face(m1,n1)) then

                number_density = rho_mol_bulk * n_avog
                radius_pol     = (3./4./pi/number_density)**(1./3.) * 1.d10

                if (n1.eq.1) then
                    r_center_surf = xc(m1,k1) - box_lo(m1) + wall_distance
                elseif (n1.eq.2) then
                    r_center_surf = box_hi(m1) - xc(m1,k1) + wall_distance
                endif

                r_surf = r_center_surf - radius_pol

                call hamaker_sphere_plate(r_surf, radius_pol, sigma_pol, sigma_plate(1), A_pol, A_plate(1), Urep, Uatt)

                Urep = Urep/(boltz_const_Joule_K*Temp)
                Uatt = Uatt/(boltz_const_Joule_K*Temp)

                Ufield(k1) = Ufield(k1) + Urep + Uatt

                if (Ufield(k1).ne.Ufield(k1)) then
                    write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". &
                                        & NaN was changed to ",E16.9)') r_surf, Ufield(k1)
                    call exit_with_error(1,2,1,ERROR_MESSAGE)
                endif

                if (r_surf.lt.0.d0) then
                    write(ERROR_MESSAGE,'("Hamaker distance smaller than zero! (",E16.9,").")') r_surf
                    call exit_with_error(1,2,1,ERROR_MESSAGE)
                endif

                write(211,'(E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3)') r_center_surf, r_surf, Uatt, Urep, Urep + Uatt
            endif
        enddo
    enddo

   !loop over all nanoparticle faces
   number_density   = rho_mol_bulk * n_avog
   radius_pol       = (3./4./pi/number_density)**(1./3.) * 1.d10

   do m1 = 1, n_nanopart_faces

       r_centers        = dsqrt((xc(1,k1)-center_np(1,m1))**2 + (xc(2,k1)-center_np(2,m1))**2 + (xc(3,k1)-center_np(3,m1))**2)
       radius_np_actual = radius_np(m1) - wall_distance
       r_surf           = r_centers - radius_pol - radius_np_actual

       call hamaker_sphere_sphere(r_surf, radius_pol, radius_np_actual, sigma_pol, sigma_np(m1), A_pol, A_np(m1), Urep, Uatt)

       Urep = Urep/(boltz_const_Joule_K*Temp)
       Uatt = Uatt/(boltz_const_Joule_K*Temp)

       Ufield(k1) = Ufield(k1) + Urep + Uatt

       if (Ufield(k1).ne.Ufield(k1)) then
           write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". &
                                & NaN was changed to ",E16.9)') r_surf, Ufield(k1)
           call exit_with_error(1,2,1,ERROR_MESSAGE)
       endif

       write(211,'(E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3)')  r_centers, r_surf, Uatt, Urep, Urep + Uatt
   enddo
enddo

close(211)

!read field
if (field_init_scheme.eq.0) then
    wa = 0.d0

elseif (field_init_scheme.eq.1) then
    write(iow,'(A40,5X,A12)')adjl("Reading field from file:",40),field_in_filename
    write(6  ,'(A40,5X,A12)')adjl("Reading field from file:",40),field_in_filename

    inquire(file = field_in_filename, exist = file_exists)

    if (file_exists) then
        open(unit=655, file = field_in_filename, Form='unformatted')
    else
        write(ERROR_MESSAGE,'("File ",A15," does not exist!")')field_in_filename
        call exit_with_error(1,1,1,ERROR_MESSAGE)
    endif

    read(655) wa
    close(655)

elseif (field_init_scheme.eq.2) then
    do k1 = 1, numnp
        if (elem_in_q0_face(k1)) then
            wa(k1) = -kapa
        else
            wa(k1) = 0.d0
        endif
    enddo
endif
!------------------------------------------------------------------------------------------------------!
end subroutine init_field
