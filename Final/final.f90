Program final
	Implicit none
	real*8 :: infty, x, dx, e, de, tol
	real*8, dimension(2) :: y
	real*8, external :: f_x, g_x
	character (len=10) :: filename
	integer :: n=1

! parameters
	tol = epsilon(0.0d0) ! tolerance of 0
	e = -8.0d0  ! lower bound energy level
	de = 1.d-2/7 ! energy resolution
	infty = 2.4d0 ! where psi(infty) is considered to be tol
	dx = de*10d0

	do while ( n.le.4 )

		x = -infty
		y(1) = tol ! initial psi at x=infty
		y(2) = tol ! initial psi' at x=infty

		write(filename , '("res",i1,".dat")') n
		open(10, file=filename)

		do while( .true. )
			call rk4(x, dx, f_x, g_x, y, e)
			if (x.lt.infty) then
				write (10,*) x, y(1), y(2)
			else
				exit
			endif
			x = x + dx
		enddo

		if (abs(y(1)).lt.tol) then
			print *, "eigenvalue", e, " psi(",x,")=",y(1)
			n = n + 1
			e = e + 0.21d-1
		endif

		e = e + de
		close(10)
	enddo

end program final


! 4th-order Runge-Kutta subroutine 
subroutine rk4(x, dx, df, dg, y, e)
	implicit none
	real*8, external :: df, dg
	real*8, intent(in) :: x, dx, e
	real*8, intent(inout), dimension(2) :: y
	real*8 :: h, k0, k1, k2, k3, l0, l1, l2, l3

	h=dx/2.0D0

	k0 = dx * df(y(1),y(2))
	l0 = dx * dg(y(1),y(2),x,e)
	k1 = dx * df(y(1)+h, y(2)+0.5d0*l0)
	l1 = dx * dg(y(1)+0.5d0*k0,y(2)+0.5d0*l0,x,e)
	k2 = dx * df(y(1)+0.5d0*k1,y(2)+0.5d0*l1)
	l2 = dx * dg(y(1)+0.5d0*k1,y(2)+0.5d0*l1,x,e)
	k3 = dx * df(y(1)+k2,y(2)+l2)
	l3 = dx * dg(y(1)+k2,y(2)+l2,x,e)
	y(1) = y(1) + (k0+2*k1+2*k2+k3)/6.0d0
	y(2) = y(2) + (l0+2*l1+2*l2+l3)/6.0d0
	Return
End subroutine rk4

! function which returns the derivatives (RHS)
real*8 function f_x(y, yp)
! dpsi/dx = psi'
	Implicit none
	real*8 ,intent(in) :: y, yp
	f_x = yp
	Return
End function f_x

real*8 function g_x(y, yp, x, e)
! d^2y/dx^2 = 2(V(x)-e)*y = (4x^4-16x^2-2e)y
	implicit none
	real*8, intent(in) :: y, yp, x, e
!	g_x = ( -4.0d0 - e) * y
!	g_x = ( 4.0d0 * x**4.0d0 - 16.0d0 * x**2.0d0 - 2.0d0*e) * y
	g_x = 2.d0 * ( (x+sqrt(2.d0))**2.d0-8 - e) * y
	Return
End function g_x