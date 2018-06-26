 	Program oscillator
 	Implicit none
c declarations
c n: number of equations, min/max in x, dist:length of x-steps
c y(1): initial position, y(2):initial velocity
	Real*8 dist, min, max, x, y(5)
	Integer	n
	n=2
        min=0.0	
        max=10.0
        dist=0.1
        y(1)=1.0
        y(2)=0.0
c open file
        Open(6, File='rk4.dat', Status='Unknown')
c do n steps of Runga-Kutta algorithm
	Do 60 x=min, max, dist
	Call rk4(x, dist, y, n)
        Write (6,*) x, y(1)
 60	Continue
c
	Close(6)
	Stop 'data saved in rk4.dat'
        End
c------------------------end of main program------------------------
c
c fourth-order Runge-Kutta subroutine 
	Subroutine rk4(x, xstep, y, n)
	Implicit none
c declarations
	Real*8 deriv, h, x, xstep, y(5) 
        Real*8 k1(5), k2(5),k3(5), k4(5), t1(5), t2(5), t3(5)
     	Integer i, n
     	h=xstep/2.0
	Do 10 i = 1,n
	   k1(i) = xstep * deriv(x, y, i)
	   t1(i) = y(i) + 0.5*k1(i)
 10	Continue
	Do 20 i = 1,n
	   k2(i) = xstep * deriv(x+h, t1, i)
	   t2(i) = y(i) + 0.5*k2(i)
 20	Continue
	Do 30 i = 1,n
	   k3(i) = xstep * deriv(x+h, t2, i)
	   t3(i) = y(i) + k3(i)
 30	Continue
	Do 40 i = 1,n
	   k4(i) = xstep * deriv(x+xstep, t3, i)
	   y(i) = y(i) + (k1(i) + (2.*(k2(i) + k3(i))) + k4(i))/6.0
 40	Continue
c 
	Return
	End
c function which returns the derivatives
	Function deriv(x, temp, i)
	Implicit none
c declarations
	Real*8 deriv, x, temp(2), omega, alpha
	Integer i
	data omega /1.d0/
	data alpha /0.5d0/
c
	If (i .EQ. 1) deriv=temp(2)
	If (i .EQ. 2) deriv=-temp(1) * omega**2 -alpha*temp(2) 
	Return
	End
