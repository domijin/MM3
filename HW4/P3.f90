Program P3
	Implicit none
	real*8 :: dt, min, max, t
	real*8, dimension(2) :: y

	! inteval
	min=0.0D0
	max=10.0D0
	
	!time step
	dt=0.1D0
	! initial position
	y(1)=1.0D0
	! initial velocity
	y(2)=0.0D0

	Open(6, File='p3.dat')
	! Runga-Kutta iteration
	Do t=min, max, dt
		Call rk4(t, dt, y, 2)
		Write (6,*) t, y(1) ! y(1): x, pos
	enddo

	Close(6)
End program P3

! 4th-order Runge-Kutta subroutine 
Subroutine rk4(t, dt, y, n)
	Implicit none
	real*8, external :: df, dg
	real*8 :: h, k0, k1, k2, k3, l0, l1, l2, l3
	real*8, intent(in) :: t, dt
	real*8, intent(inout), dimension(2) :: y
	integer, intent(in) :: n

	h=dt/2.0D0
	k0 = dt * df(y(1),y(n))
	l0 = dt * dg(y(1),y(n))
	k1 = dt * df(y(1)+h, y(n)+0.5d0*l0)
	l1 = dt * dg(y(1)+0.5d0*k0,y(n)+0.5d0*l0)
	k2 = dt * df(y(1)+0.5d0*k1,y(n)+0.5d0*l1)
	l2 = dt * dg(y(1)+0.5d0*k1,y(n)+0.5d0*l1)
	k3 = dt * df(y(1)+k2,y(n)+l2)
	l3 = dt * dg(y(1)+k2,y(n)+l2)
	y(1) = y(1) + (k0+2*k1+2*k2+k3)/6.0d0
	y(n) = y(n) + (l0+2*l1+2*l2+l3)/6.0d0
	Return
End subroutine rk4

! function which returns the derivatives (RHS)
real*8 function df(y,z)
! dx/dt = v(t)
	Implicit none
	real*8 ,intent(in) :: y, z
	df=z
	Return
End function df

real*8 Function dg(y,z)
! dv/dt = -w**2*x(t) -a*v
! the second term is damping term
	implicit none
	real*8, intent(in) :: y, z
	real*8 :: omega, alpha
	data omega /3.0d0/
	data alpha /0.5d0/

	dg= - omega**2.0d0 * y - alpha * z
	Return
End function dg