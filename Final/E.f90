Program final
	Implicit none
	real*8 :: infty, x, dx, e1, e2, e
	real*8, dimension(2) :: y
	real*8, external :: f_x, g_x
	character (len=10) :: filename
	integer :: n

! parameters
	infty = 4.d0 ! where psi(infty) is considered to be tol

	dx = 0.1d-2/7.d0

	! e1 = -5.3228961107477559
	! e2 = -5.3228961107462656
	! n = 1

	! e1 = -5.3034634104460840
	! e2 = -5.3034634104460832
	! n = 2


	e = -5.3034634104460840
	n = 2

	x = -infty

	y(1) = epsilon(0.d0) ! initial psi at x=infty
	y(2) = epsilon(0.d0) ! initial psi' at x=infty


	write(filename , '("res",i1,".dat")') n
	open(10, file=filename,action="write",status="replace")

	do while( x.le.infty )
		call rk4(x, dx, f_x, g_x, y, e)
		write (10,*) x, y(1)
		x = x + dx
	enddo

	close(10)
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
	g_x = ( 4.0d0 * x**4.0d0 - 16.0d0 * x**2.0d0 - 2.0d0*e) * y
!	g_x = 2.d0 * ( (x+sqrt(2.d0))**2.d0-8 - e) * y
	Return
End function g_x