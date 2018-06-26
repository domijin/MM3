Program proj2
	Implicit none
	real*8 :: t, dt, v, theta0, theta, pi, offset, alpha
	real*8, dimension(2) :: x, y
	real*8, external :: f_x, g_x, g_y
	integer :: i
	character (len=10) :: filename

	pi = 4.0d0*atan(1.0d0)

	! read setting
	print *, "Initial angle in degree:"
	read *, theta0

	do i=1,2
		write(filename , '("res",i1,".dat")') i
		theta = theta0 /180.0d0 * pi
		alpha = 5.0d-2 / 30.0d0 * (i-1.d0) ! k/m

		! initial conditions
		v = 100.0D0
		t = 0.0d0
		dt = 0.02D0

		x(1) = 0.0D0
		x(2) = v * cos(theta)

		y(1) = EPSILON(0.0d0)
		y(2) = v * sin(theta)

		open(10, file=filename)

		! Runga-Kutta iteration
		do while( .true. )
			t = t + dt
			call rk4(dt, f_x, g_x, x, v, theta, alpha)
			call rk4(dt, f_x, g_y, y, v, theta, alpha)
			v = sqrt(x(2)**2 + y(2)**2)
			theta = atan(y(2)/x(2))
			if (y(1).gt.0) then
				write (10,*) t, x(1), y(1), x(2), y(2), v, theta / pi * 180 
			else
				exit
			endif
		enddo

		offset = y(1) / y(2)
		t = t - offset
		x(1) = x(1) - offset * x(2)
		y(1) = y(1) - offset * y(2) 
		write (10,*) t, x(1), y(1), x(2), y(2), v, theta / pi * 180
		close(10)
		print *, 'flying time:', t, 's; distance:', x(1), 'meter'
	enddo
End program proj2

! 4th-order Runge-Kutta subroutine 
subroutine rk4(dt, df, dg, y, v, theta, alpha)
	implicit none
	real*8, external :: df, dg
	real*8, intent(in) :: dt, alpha
	real*8, intent(inout) :: v, theta
	real*8, intent(inout), dimension(2) :: y
	real*8 :: h, k0, k1, k2, k3, l0, l1, l2, l3

	h=dt/2.0D0

	k0 = dt * df(y(1),y(2))
	l0 = dt * dg(y(1),y(2),v,theta,alpha)
	k1 = dt * df(y(1)+h, y(2)+0.5d0*l0)
	l1 = dt * dg(y(1)+0.5d0*k0,y(2)+0.5d0*l0,v,theta,alpha)
	k2 = dt * df(y(1)+0.5d0*k1,y(2)+0.5d0*l1)
	l2 = dt * dg(y(1)+0.5d0*k1,y(2)+0.5d0*l1,v,theta,alpha)
	k3 = dt * df(y(1)+k2,y(2)+l2)
	l3 = dt * dg(y(1)+k2,y(2)+l2,v,theta,alpha)
	y(1) = y(1) + (k0+2*k1+2*k2+k3)/6.0d0
	y(2) = y(2) + (l0+2*l1+2*l2+l3)/6.0d0
	Return
End subroutine rk4

! function which returns the derivatives (RHS)
real*8 function f_x(a, b)
! dx/dt = v(t)
	Implicit none
	real*8 ,intent(in) :: a, b
	f_x = b
	Return
End function f_x

real*8 function g_x(a, b, v, theta, alpha)
	implicit none
	real*8, intent(in) :: a, b, v, theta, alpha
	g_x = - alpha * v**2 * cos(theta)
!	g_x = 0
	Return
End function g_x

real*8 function g_y(a, b, v, theta, alpha)
	implicit none
	real*8, intent(in) :: a, b, v, theta, alpha
	real*8 :: g
	g = 9.8d0
	g_y = -g - alpha * v**2 * sin(theta)
!	g_y = -g
	Return
End function g_y