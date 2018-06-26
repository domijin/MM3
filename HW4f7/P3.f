 	Program P3
 	Implicit none
	Real*8 dt, min, max, t, y(5)
c n: number of equations (the 'dimension' of the standard form)
C starting and end points
        min=0.0D0
        max=10.0D0
C time step
        dt=0.1D0
C initial position
        y(1)=1.0D0
C initial velocity
        y(2)=0.0D0
c open file to output
        Open(6, File='p3.dat')
c do n steps of Runga-Kutta algorithm
	Do t=min, max, dt
	  Call rk4(t, dt, y, 2)
          Write (6,*) t, y(1)
	enddo
c
	Close(6)
	Stop 
        End
c------------------------end of main program------------------------
c
	Subroutine rk4(t, dt, y, n)
	Implicit none
	Real*8 df, dg, h, t, dt, y(5)
        Real*8 k0, k1, k2, k3
        Real*8 l0, l1, l2, l3
        integer n
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
	End

c function which returns the derivatives (RHS)
	Function df(y,z)
C 2 function components: dx/dt = v(t),   dv/dt = -w**2*x(t) -a*v
C the second term is damping term
	Implicit none
c declarations
	Real*8 df, y, z
	df=z
	Return
	End

	Function dg(y,z)
C 2 function components: dx/dt = v(t),   dv/dt = -w**2*x(t) -a*v
C the second term is damping term
	implicit none
c declarations
	Real*8 dg, y, z, omega, alpha
	data omega /3.0d0/
	data alpha /0.5d0/
c
	dg= - omega**2.0d0 * y - alpha * z
	Return
	End