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
real(8)                                :: number_density = 0.d0, distance = 0.d0, radius_pol = 0.d0
real(8)                                :: Urep = 0.d0, Uatt = 0.d0
!------------------------------------------------------------------------------------------------------!
open(unit=211, file = usolid)

Ufield = 0.d0

do k1 = 1, numnp
    !loop over all dirichlet faces
    do m1 = 1, 3
        do n1 = 1, 2
            if (is_dir_face(m1,n1)) then

                number_density = rho_0 * avogadro_constant
                radius_pol     = (3./4./pi/number_density)**(1./3.)
                wall_distance  = radius_pol

                if (n1.eq.1) then
                    distance = xc(m1,k1) - box_lo(m1) - radius_pol + wall_distance
                elseif (n1.eq.2) then
                    distance = box_hi(m1) - xc(m1,k1) - radius_pol + wall_distance
                endif

                call hamaker_sphere_plate(distance, radius_pol, sigma_pol, sigma_plate(1), A_pol, A_plate(1), Urep, Uatt)

                Urep = Urep/(boltz_const_Joule_K*Temp)
                Uatt = Uatt/(boltz_const_Joule_K*Temp)

                if ((Urep+Uatt).gt.0.) then
                    Urep = 0.d0
                    Uatt = 0.d0
                endif

                Ufield(k1) = Ufield(k1) + Urep + Uatt

                if (Ufield(k1).ne.Ufield(k1)) then
                    write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". &
                                        & NaN was changed to ",E16.9)') distance, Ufield(k1)
                    call exit_with_error(1,2,1,ERROR_MESSAGE)
                endif

                if (distance.lt.0.d0) then
                    write(ERROR_MESSAGE,'("Hamaker distance smaller than zero! (",E16.9,").")') distance
                    call exit_with_error(1,2,1,ERROR_MESSAGE)
                endif

                write(211,'(E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3)') distance, Uatt, Urep, Urep + Uatt
            endif
        enddo
    enddo

   !loop over all nanoparticle faces
   do m1 = 1, n_nanopart_faces

       number_density = rho_0 * avogadro_constant
       radius_pol     = (3./4./pi/number_density)**(1./3.)
       radius_np(m1)  = radius_np(m1)/1.e+10
       distance       = dsqrt((xc(1,k1)-center_np(1,m1))**2 + (xc(2,k1)-center_np(2,m1))**2 + (xc(3,k1)-center_np(3,m1))**2) &
                   &  - radius_pol - radius_np(m1)

       call hamaker_sphere_sphere(distance, radius_pol, radius_np(m1), sigma_pol, sigma_np(m1), A_pol, A_np(m1), Urep, Uatt)

       Urep = Urep/(boltz_const_Joule_K*Temp)
       Uatt = Uatt/(boltz_const_Joule_K*Temp)

       if ((Urep+Uatt).gt.0.) then
           Urep = 0.d0
           Uatt = 0.d0
       endif

       Ufield(k1) = Ufield(k1) + Urep + Uatt

       if (Ufield(k1).ne.Ufield(k1)) then
           write(ERROR_MESSAGE,'("Hamaker assumed a NaN value for x = ",E16.9,". &
                                & NaN was changed to ",E16.9)') distance, Ufield(k1)
           call exit_with_error(1,2,1,ERROR_MESSAGE)
       endif

       write(211,'(E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3)') distance, Uatt, Urep, Urep + Uatt
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
