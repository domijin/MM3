Program P1
	Implicit none
	real*8 :: dt, min, max, t
	real*8, dimension(2) :: y
	integer :: n
	parameter (n=2) ! y 1: pos; 2: vel; 

	! inteval
	min=0.0D0
	max=200.0D0
	! time step
	dt=0.1D0
	! initial position
	y(1)=0.0D0
	! initial velocity
	y(2)=4.0D0

	! print the exact result
	print *, sqrt(4.0D0**2.D0+2.d0*400*max/70)
	Open(6, File='P1.dat')

	! Runga-Kutta iteration
	Do t=min, max, dt
		Call rk2(t, dt, y, n)
		Write (6,*) t, y(2) ! y(2): vel
	enddo

	Close(6)
End program P1

! second-order Runge-Kutta subroutine 
! will update y(1) & y(2) based on n
subroutine rk2(t, dt, y, n)
	Implicit none
	real*8, external :: deriv
	real*8, intent(in) :: t, dt
	real*8, intent(inout), dimension(2) :: y
	real*8 :: h
	real*8, dimension(2) :: k1, k2, t1
	integer :: i,n

	h=dt/2.0D0
	Do i = 1,n
	   k1(i) = dt * deriv(y, i)
	   t1(i) = y(i) + 0.5D0*k1(i)
	enddo
	Do i = 1,n
	   k2(i) = dt * deriv(t1, i)
	   y(i) = y(i) + k2(i)
	enddo
	Return
End subroutine rk2

! function which returns the derivatives (RHS)
real*8 function deriv(temp, i)
	Implicit none
! declarations
	Real*8 :: omega
	real*8, dimension(2), intent(inout):: temp
	Integer :: i
	data omega /5.7143d0/
	!   dx/dt=v
	!   dv/dt=P/(mv)= 400/70/v=5.7143/v
	If (i .EQ. 1) deriv=temp(2)
	If (i .EQ. 2) deriv=omega / temp(2)
	Return
End function deriv
