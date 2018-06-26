 	Program oscillator
 	Implicit none
c declarations
	Real*8 dt, min, max, t, y(5)
c n: number of equations (the 'dimension' of the standard form)
	Integer	n
C number of component
        parameter (n=2)
C starting and end points
        min=0.0D0	
        max=10.0D0
C time step
        dt=0.1D0
C initial position
        y(1)=0.2D0
C initial velocity
        y(2)=0.0D0

c open file to output
        Open(6, File='rk2.dat')
c do n steps of Runga-Kutta algorithm
	Do t=min, max, dt
	  Call rk2(t, dt, y, n)
          Write (6,*) t, y(1)
   	enddo
c
	Close(6)
	Stop 'data saved in rk2.dat'
        End
c------------------------end of main program------------------------
c
c second-order Runge-Kutta subroutine 
	Subroutine rk2(t, dt, y, n)
	Implicit none
c declarations
	Real*8 deriv, h, t, dt, y(1) 
        Real*8 k1(5), k2(5), t1(5)
     	Integer i, n
     	h=dt/2.0D0
	Do i = 1,n
	   k1(i) = dt * deriv(t, y, i)
	   t1(i) = y(i) + 0.5D0*k1(i)
   	enddo
	Do i = 1,n
	   k2(i) = dt * deriv(t+h, t1, i)
	   y(i) = y(i) + k2(i)
   	enddo
c 
	Return
	End

c function which returns the derivatives
	Function deriv(x, temp, i)
C 2 function components: dx/dt = v(t),   dv/dt = -w**2*x(t) -a*v
C the second term is damping term
	Implicit none
c declarations
	Real*8 deriv, x, temp(2), omega, alpha
	Integer i
	data omega /3.13d0/
	data alpha /0.5d0/
c
	If (i .EQ. 1) deriv=temp(2)
	If (i .EQ. 2) deriv=-temp(1) * omega**2
	Return
	End
