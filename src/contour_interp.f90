!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine interp_lagrange(m, data_num, t_data, p_data, interp_num, t_interp, p_interp)
implicit none

integer :: m, data_num, interp_num

real(8) :: l_interp(data_num,interp_num)
real(8) :: p_data(m,data_num)
real(8) :: p_interp(m,interp_num)
real(8) :: t_data(data_num)
real(8) :: t_interp(interp_num)

!Evaluate the DATA_NUM Lagrange polynomials associated with T_DATA(1:DATA_NUM)
!for the interpolation points T_INTERP(1:INTERP_NUM).
call lagrange_value(data_num, t_data, interp_num, t_interp, l_interp)
!Multiply P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
!to get P_INTERP(1:M,1:INTERP_NUM).
p_interp(1:m,1:interp_num) = MATMUL (p_data(1:m,1:data_num), l_interp(1:data_num,1:interp_num))
return
end subroutine interp_lagrange


subroutine interp_linear (m, data_num, t_data, p_data, interp_num, t_interp, p_interp)
implicit none

integer :: m, data_num, interp_num
integer :: interp, left, right

logical :: r8vec_ascends_strictly

real(8) :: p_data(m,data_num)
real(8) :: p_interp(m,interp_num)
real(8) :: t
real(8) :: t_data(data_num)
real(8) :: t_interp(interp_num)

if (.not.r8vec_ascends_strictly(data_num, t_data)) then
  write(*,'(A)') ' '
  write(*,'(A)') "INTERP_LINEAR - Fatal error!"
  write(*,'(A)') "  Independent variable array T_DATA is not strictly increasing."
  STOP 1
endif

do interp = 1, interp_num
    t = t_interp(interp)
    !Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
    !nearest to, TVAL.
    call r8vec_bracket(data_num, t_data, t, left, right)

    p_interp(1:m,interp) = ((t_data(right)-t) * p_data(1:m,left) + (t-t_data(left)) * p_data(1:m,right)) &
                                                               & / (t_data(right)-t_data(left))
enddo
return
end subroutine interp_linear


subroutine lagrange_value(data_num, t_data, interp_num, t_interp, l_interp)
implicit none

integer :: data_num, interp_num
integer :: i, j

real(8) :: l_interp(data_num,interp_num)
real(8) :: t_data(data_num)
real(8) :: t_interp(interp_num)

!Evaluate the polynomial.
l_interp(1:data_num,1:interp_num) = 1.0D+00

do i = 1, data_num
    do j = 1, data_num
        if (j/=i) then
            l_interp(i,1:interp_num) = l_interp(i,1:interp_num) &
                                       * ( t_interp(1:interp_num) - t_data(j) ) / ( t_data(i) - t_data(j) )
        endif
    enddo
enddo
return
end subroutine lagrange_value


logical function r8vec_ascends_strictly(n, x)
implicit none

integer :: n
integer :: i

real(8) :: x(n)

do i = 1, n - 1
    if (x(i+1)<=x(i)) then
        r8vec_ascends_strictly = .false.
        return
    endif
enddo
r8vec_ascends_strictly = .true.
return
end function r8vec_ascends_strictly


subroutine r8vec_bracket(n, x, xval, left, right)
implicit none

integer :: n
integer :: left, right
integer :: i

real(8) :: x(n)
real(8) :: xval

do i = 2, n - 1
    if (xval<x(i)) then
        left  = i - 1
        right = i
        return
    endif
enddo
left  = n - 1
right = n
return
end subroutine r8vec_bracket


subroutine timestamp( )
implicit none

character(len=8)                           :: ampm
character(len=9), parameter, dimension(12) :: month = (/ &
  "January  ", "February ", "March    ", "April    ", &
  "May      ", "June     ", "July     ", "August   ", &
  "September", "October  ", "November ", "December " /)

integer :: d, h, m, mm, n, s, y
integer :: values(8)

call DATE_AND_TIME(values=values)

y  = values(1)
m  = values(2)
d  = values(3)
h  = values(5)
n  = values(6)
s  = values(7)
mm = values(8)

if (h<12) then
    ampm = "AM"
elseif (h==12) then
    if (n==0.and.s==0) then
        ampm = "Noon"
    else
        ampm = "PM"
    endif
else
    h = h - 12
    if (h<12) then
        ampm = "PM"
    else if ( h == 12 ) then
        if (n==0.and.s==0) then
            ampm = "Midnight"
        else
            ampm = "AM"
        endif
    endif
endif
write(*,'(I2,1X,A,1X,I4,2X,I2,A1,I2.2,A1,I2.2,A1,I3.3,1X,A)') d, TRIM(month(m)), y, h, ':', n, ':', s, '.', mm, TRIM(ampm)
return
end subroutine timestamp
