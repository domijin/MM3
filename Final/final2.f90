Program final2
	Implicit none
	real*8 :: infty, x, dx, e1, e2, e, tol, phi1, phi2
	real*8, dimension(2) :: y1, y2, y
	real*8, external :: f_x, g_x
	character (len=10) :: filename
	integer :: n=1

! parameters
	tol = epsilon(0.0d0) ! tolerance of 0
	infty = 2.4d0 ! where psi(infty) is considered to be tol
	e1 = -8.0d0  ! lower bound energy level

	do while ( n.le.4 ) ! get first 4 lowest energy 

		x = -infty
		dx = infty/1.0d4
		y1(1) = tol
		y1(2) = tol
		do while( x.le.infty ) ! integrate to get psi1(infty)
			call rk4(x, dx, f_x, g_x, y1, e1)
			x = x + dx
		enddo
		phi1 = y1(1)
		phi2 = y1(1)
		e2 = e1

		print *, "e1 and phi1", e1, phi1

		do while ( phi1*phi2.gt.0 ) ! find e2
			x = -infty
			e2 = e2 + 3.d0/7.d0
			dx = infty/1.0d5
			y2(1) = tol
			y2(2) = tol

			do while( x.le.infty ) ! integrate to get psi2(infty)
				call rk4(x, dx, f_x, g_x, y2, e2)
				x = x + dx
			enddo
			phi2 = y2(1)

			print *, "finding e2", e2, phi2
		enddo

		y(1) = 1.0
		do while (abs(phi1-phi2).gt.tol) ! find e by updating e1 & e2
			x = -infty
			e = (e1 + e2)/2.d0
			dx = ( e2 - e1 )/9.d0
			y1(1) = tol
			y2(1) = tol
			y(1) = tol
			y1(2) = tol
			y2(2) = tol
			y(2) = tol

			write(filename , '("res",i1,".dat")') n
			open(10, file=filename,action="write",status="replace")
			do while( x.le.infty )

				call rk4(x, dx, f_x, g_x, y1, e1)
				call rk4(x, dx, f_x, g_x, y2, e2)
				call rk4(x, dx, f_x, g_x, y, e)
				write (10,*) x, y(1), y(2)
				x = x + dx
			enddo
			phi1 = y1(1)
			phi2 = y2(1)

			if (phi1*y(1).lt.0) then
				e2 = e
			else
				e1 = e
			endif

			close(10)

!			print *, "updating e1, e2", e1, e2, y(1)
		enddo
		e1 = e + 0.1d0/7.d0

		print *, "eigenvalue",n,"is:", e
		n = n + 1
	enddo

end program final2


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
	Return
End function g_x