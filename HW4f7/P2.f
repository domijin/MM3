 	Program P2
 	Implicit none
	Real*8 dt, min, max, t, y(5)
c n: number of equations (the 'dimension' of the standard form)
	Integer	n
C number of component
        parameter (n=2)
C 	 y 1: pos; 2: vel; 
C starting and end points
        min=0.0D0
        max=200.0D0
C time step
        dt=0.1D0
C initial position
        y(1)=0.0D0
C initial velocity
        y(2)=4.0D0
c open file to output
        Open(6, File='p2.dat')
c do n steps of Runga-Kutta algorithm
	Do t=min, max, dt
	  Call rk2(t, dt, y, n)
          Write (6,*) t, y(2)
   	enddo
c
	Close(6)
	Stop 
        End
c------------------------end of main program------------------------
c
c second-order Runge-Kutta subroutine 
c will update y(1) & y(2) based on n=2
	Subroutine rk2(t, dt, y, n)
	Implicit none
	Real*8 deriv, h, t, dt, y(1) 
        Real*8 k1(5), k2(5), t1(5)
     	Integer i, n

     	h=dt/2.0D0
	Do i = 1,n
	   k1(i) = dt * deriv(y, i)
	   t1(i) = y(i) + 0.5D0*k1(i)
   	enddo
	Do i = 1,n
	   k2(i) = dt * deriv(t1, i)
	   y(i) = y(i) + k2(i)
   	enddo
c 
	Return
	End

c function which returns the derivatives (RHS)
	Function deriv( temp, i)
C 2 function components: dx/dt = v(t),   dv/dt = -w**2*x(t) -a*v
C the second term is damping term
c  i=1, return temp(2)->y(2): vel
c  i=2, return -w**2*temp(1)-> -w^2*y(1): pos
	Implicit none
c declarations
	Real*8 deriv, temp(2), omega, alpha
	Integer i
	data omega /5.7143d0/
c   dv/dt=P/(mv)= 400/70/v=5.7143/v
c   dx/dt=v
	data alpha /0.15d0/
c
	If (i .EQ. 1) deriv=temp(2)
	If (i .EQ. 2) deriv=omega / temp(2) - alpha* temp(2)**2
	Return
	End
