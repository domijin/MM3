 	Program oscillator
 	Implicit none
c declarations
c n: number of equations
	Real*8 dist, min, max, t, y(5)
	Integer	n
	n=2
C starting and end points
        min=0.0D0	
        max=10.0D0
C time step
        dist=0.1D0
C initial position
        y(1)=0.2
C initial velocity
        y(2)=0.0

c open file to output
        Open(6, File='rk4.dat')
c do n steps of Runga-Kutta algorithm
	Do t=min, max, dist
	  Call rk4(t, dist, y, n)
          Write (6,100) t, y(1)
 	enddo
c
	Close(6)
	Stop 'data saved in rk4.dat'
100     format(2f14.6)
        End
c------------------------end of main program------------------------
c
c fourth-order Runge-Kutta subroutine 
	Subroutine rk4(x, xstep, y, n)
	Implicit none
c declarations
	Real*8 func, h, x, xstep, y(1) 
        Real*8 k1(5), k2(5),k3(5), k4(5), t1(5), t2(5), t3(5)
     	Integer i, n
C half step size
     	h=xstep/2.D0
	Do i = 1,n
	   k1(i) = xstep * func(x, y, i)
	   t1(i) = y(i) + 0.5*k1(i)
   	enddo
	Do i = 1,n
	   k2(i) = xstep * func(x+h, t1, i)
	   t2(i) = y(i) + 0.5*k2(i)
   	enddo
	Do i = 1,n
	   k3(i) = xstep * func(x+h, t2, i)
	   t3(i) = y(i) + k3(i)
   	enddo
	Do i = 1,n
	   k4(i) = xstep * func(x+xstep, t3, i)
	   y(i) = y(i) + (k1(i) + (2.*(k2(i) + k3(i))) + k4(i))/6.0
   	enddo
c 
	Return
	End

c function which returns the func
	Function func(t, temp, i)
C 2 function components: dx/dt = v(t),   dv/dt = -w**2*x(t) -a*v
C the second term is damping term
	Implicit none
c declarations
	Real*8 func, t, temp(2), omega, alpha
	Integer i
	data omega /3.13d0/
	data alpha /0.5d0/
c
	If (i .EQ. 1) func=temp(2)
	If (i .EQ. 2) func=-temp(1) * omega**2
	Return
	End
