subroutine init_field(field_in_filename, Ufield, wa)
!------------------------------------------------------------------------------------------------------!
use parser_vars
use geometry
use error_handing
use write_helper
use force_fields
use iofiles
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer :: k1, m1, n1

real(8), intent(out), dimension(numnp) :: Ufield, wa
real(8)                                :: distance, r12, Urep, Uatt, Utot
real(8)                                :: sigma_plate_temp, A_plate_temp
!------------------------------------------------------------------------------------------------------!
open(unit=211, file = usolid)

Ufield = 0.d0

do k1 = 1, numnp
    !loop over all dirichlet faces
    do m1 = 1, 3
        do n1 = 1, 2
            if (is_dir_face(m1,n1)) then

                sigma_plate_temp        = sigma_plate(1)
                A_plate_temp            = A_plate(1)

                if (n1.eq.1) then
                    distance = xc(m1,k1) - box_lo(m1)
                elseif (n1.eq.2) then
                    distance = box_hi(m1) - xc(m1,k1)
                endif

                call hamaker_sphere_plate(Temp, distance, rho_0, sigma_pol, sigma_plate_temp, A_pol, &
    &                                     A_plate_temp, r12, Urep, Uatt, Utot)
                Ufield(k1) = Ufield(k1) + Utot

                if (Ufield(k1).ne.Ufield(k1)) then
                    write(ERROR_MESSAGE,'(''Hamaker assumed a NaN value for x = '',E16.9,''. &
                                        & NaN was changed to '',E16.9)') distance, Ufield(k1)
                    call exit_with_error(1,2,1,ERROR_MESSAGE)
                endif

                if (distance.lt.0.d0) then
                    write(ERROR_MESSAGE,'(''Hamaker distance smaller than zero! ('',E16.9,'').'')') distance
                    call exit_with_error(1,2,1,ERROR_MESSAGE)
                endif

                write(211,'(E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3)') distance, r12, Uatt, Urep, Utot
            endif
        enddo
    enddo

   !loop over all nanoparticle faces
   do m1 = 1, n_nanopart_faces

           distance = dsqrt( (xc(1,k1)-center_np(1,m1))**2+(xc(2,k1)-center_np(2,m1))**2 &
   &                        +(xc(3,k1)-center_np(3,m1))**2 )

           call hamaker_sphere_sphere(Temp, distance, radius_np(m1), rho_0, sigma_pol, sigma_np(m1), A_pol, A_np(m1), &
   &                                  r12, Urep, Uatt, Utot)
           Ufield(k1) = Ufield(k1) + Utot

           if (Ufield(k1).ne.Ufield(k1)) then
               write(ERROR_MESSAGE,'(''Hamaker assumed a NaN value for x = '',E16.9,''. &
                                   & NaN was changed to '',E16.9)') distance, Ufield(k1)
               call exit_with_error(1,2,1,ERROR_MESSAGE)
           endif

           write(211,'(E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3,2X,E17.9E3)') distance, r12, Uatt, Urep, Utot
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
