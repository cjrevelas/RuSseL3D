subroutine interp_lagrange(m, data_num, t_data, p_data, interp_num, t_interp, p_interp)
implicit none

integer ( kind = 4 ) data_num
integer ( kind = 4 ) m
integer ( kind = 4 ) interp_num

real ( kind = 8 ) l_interp(data_num,interp_num)
real ( kind = 8 ) p_data(m,data_num)
real ( kind = 8 ) p_interp(m,interp_num)
real ( kind = 8 ) t_data(data_num)
real ( kind = 8 ) t_interp(interp_num)

!Evaluate the DATA_NUM Lagrange polynomials associated with T_DATA(1:DATA_NUM)
!for the interpolation points T_INTERP(1:INTERP_NUM).

call lagrange_value(data_num, t_data, interp_num, t_interp, l_interp)

!Multiply P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
!to get P_INTERP(1:M,1:INTERP_NUM).

p_interp(1:m,1:interp_num) = matmul (p_data(1:m,1:data_num), l_interp(1:data_num,1:interp_num))

return
end subroutine interp_lagrange





subroutine interp_linear (m, data_num, t_data, p_data, interp_num, t_interp, p_interp)
implicit none

integer ( kind = 4 ) data_num
integer ( kind = 4 ) m
integer ( kind = 4 ) interp_num

integer ( kind = 4 ) interp
integer ( kind = 4 ) left
real ( kind = 8 ) p_data(m,data_num)
real ( kind = 8 ) p_interp(m,interp_num)
logical r8vec_ascends_strictly
integer ( kind = 4 ) right
real ( kind = 8 ) t
real ( kind = 8 ) t_data(data_num)
real ( kind = 8 ) t_interp(interp_num)



if (.not.r8vec_ascends_strictly(data_num, t_data)) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTERP_LINEAR - Fatal error!'
  write ( *, '(a)' ) &
    '  Independent variable array T_DATA is not strictly increasing.'
  stop 1
end if

do interp = 1, interp_num

  t = t_interp(interp)

!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.

  call r8vec_bracket(data_num, t_data, t, left, right)

  p_interp(1:m,interp) = ((t_data(right)-t) * p_data(1:m,left) + (t-t_data(left)) * p_data(1:m,right)) &
                                                                        & / (t_data(right)-t_data(left))

end do

return
end subroutine interp_linear





subroutine lagrange_value(data_num, t_data, interp_num, t_interp, l_interp)
implicit none

integer ( kind = 4 ) data_num
integer ( kind = 4 ) interp_num

integer ( kind = 4 ) i
integer ( kind = 4 ) j
real ( kind = 8 ) l_interp(data_num,interp_num)
real ( kind = 8 ) t_data(data_num)
real ( kind = 8 ) t_interp(interp_num)

!  Evaluate the polynomial.

l_interp(1:data_num,1:interp_num) = 1.0D+00

do i = 1, data_num

  do j = 1, data_num

    if ( j /= i ) then

      l_interp(i,1:interp_num) = l_interp(i,1:interp_num) &
        * ( t_interp(1:interp_num) - t_data(j) ) / ( t_data(i) - t_data(j) )

    end if

  end do

end do

return
end subroutine lagrange_value





function r8vec_ascends_strictly(n, x)
implicit none

integer ( kind = 4 ) n

integer ( kind = 4 ) i
logical r8vec_ascends_strictly
real ( kind = 8 ) x(n)

do i = 1, n - 1
  if ( x(i+1) <= x(i) ) then
    r8vec_ascends_strictly = .false.
    return
  end if
end do

r8vec_ascends_strictly = .true.

return
end function r8vec_ascends_strictly





subroutine r8vec_bracket(n, x, xval, left, right)
implicit none

integer ( kind = 4 ) n

integer ( kind = 4 ) i
integer ( kind = 4 ) left
integer ( kind = 4 ) right
real ( kind = 8 ) x(n)
real ( kind = 8 ) xval

do i = 2, n - 1

  if ( xval < x(i) ) then
    left = i - 1
    right = i
    return
  end if

end do

left = n - 1
right = n

return
end subroutine r8vec_bracket





subroutine timestamp( )
implicit none

character ( len = 8 ) ampm
integer ( kind = 4 ) d
integer ( kind = 4 ) h
integer ( kind = 4 ) m
integer ( kind = 4 ) mm
character ( len = 9 ), parameter, dimension(12) :: month = (/ &
  'January  ', 'February ', 'March    ', 'April    ', &
  'May      ', 'June     ', 'July     ', 'August   ', &
  'September', 'October  ', 'November ', 'December ' /)
integer ( kind = 4 ) n
integer ( kind = 4 ) s
integer ( kind = 4 ) values(8)
integer ( kind = 4 ) y

call date_and_time ( values = values )

y = values(1)
m = values(2)
d = values(3)
h = values(5)
n = values(6)
s = values(7)
mm = values(8)

if ( h < 12 ) then
  ampm = 'AM'
else if ( h == 12 ) then
  if ( n == 0 .and. s == 0 ) then
    ampm = 'Noon'
  else
    ampm = 'PM'
  end if
else
  h = h - 12
  if ( h < 12 ) then
    ampm = 'PM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Midnight'
    else
      ampm = 'AM'
    end if
  end if
end if

write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
  d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

return
end subroutine timestamp
