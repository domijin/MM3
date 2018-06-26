# HW4

Path to source code: `/home/d/dx/dxj4360/HW4`
## P1: differential equation by rk2
Differential equation: $\frac{dv}{dt}=\frac{P}{mv}$, integrate from $[0,200]$ with $v_0=4.0m/s$, $P=400Watt$, $m=70kg$

The exact solution is: $v=\sqrt{v_0^2+2Pt/m}$

### modification of the DE function
$\frac{dv}{dt}=\frac{P}{mv}=\omega/v$

$\omega=\frac{P}{m}=5.7143$

```fortran
Function deriv(x, temp, i)
	Implicit none
	Real*8 deriv, x, temp(2), omega, alpha
	Integer i
	data omega /5.7143d0/
	If (i .EQ. 1) deriv=temp(2)
	If (i .EQ. 2) deriv=omega / temp(2)
Return
```

###Plot of the Result
![](https://s31.postimg.org/pgl225ivf/Screenshot_2016_07_18_17_32_00.png)

The final result is

```
   0.0000000000000000        4.1403512219144201
  0.10000000000000001        4.2761034937351496
  0.20000000000000001        4.4076808606527766
  0.30000000000000004        4.5354460088889015
  0.40000000000000002        4.6597120097890157
  0.50000000000000000        4.7807513250099660
  0.59999999999999998        4.8988028147021838
  0.69999999999999996        5.0140772643314993
  0.79999999999999993        5.1267617955466305
  0.89999999999999991        5.2370234245831986
  ...
   199.09999999999297        47.880809961041173
   199.19999999999297        47.892742900230651
   199.29999999999296        47.904672866953852
   199.39999999999296        47.916599863430982
   199.49999999999295        47.928523891879465
   199.59999999999295        47.940444954513978
   199.69999999999294        47.952363053546449
   199.79999999999293        47.964278191186054
   199.89999999999293        47.976190369639220
   199.99999999999292        47.988099591109645
```

Which is pretty close to the exact result`47.976184568119692`.

## P2: with air-grad term
### modification of the DE function
$\frac{dv}{dt}=\frac{P}{mv}-Av^2=\omega/v-\alpha v^2$

$\omega=\frac{P}{m}=5.7143$

$\alpha=-0.15$

```fortran
Function deriv(x, temp, i)
	Implicit none
	Real*8 deriv, x, temp(2), omega, alpha
	Integer i
	data omega /5.7143d0/
	data alpha /0.15d0/
	If (i .EQ. 1) deriv=temp(2)
	If (i .EQ. 2) deriv=omega / temp(2) - alpha * temp(2)**2
Return
```

###Plot of the Result, compared to P1
![](https://s32.postimg.org/3npl9pp8l/Screenshot_2016_07_18_17_45_26.png)

The final result is

```
   0.0000000000000000        3.9104066773002248
  0.10000000000000001        3.8335965544074022
  0.20000000000000001        3.7676949971424962
  0.30000000000000004        3.7111194936880341
  0.40000000000000002        3.6625285040563833
  0.50000000000000000        3.6207809805232194
  0.59999999999999998        3.5849038738538459
  0.69999999999999996        3.5540657176454782
  0.79999999999999993        3.5275549167470794
  0.89999999999999991        3.5047617352693576
  ...
   199.09999999999297        3.3647845354623893
   199.19999999999297        3.3647845354623893
   199.29999999999296        3.3647845354623893
   199.39999999999296        3.3647845354623893
   199.49999999999295        3.3647845354623893
   199.59999999999295        3.3647845354623893
   199.69999999999294        3.3647845354623893
   199.79999999999293        3.3647845354623893
   199.89999999999293        3.3647845354623893
   199.99999999999292        3.3647845354623893
```

Dramatic difference of the drag term.

## P3: 4th Runge-Kutta for damped oscillation problem
### modification of the DE function
$\frac{d^2x}{dt^2}=-\omega^2x-\alpha \frac{dx}{dt}=-\omega^2 x-\alpha v$

$\omega=3.0$

$\alpha=0.5$

```fortran
Function deriv(x, temp, i)
	Implicit none
	Real*8 deriv, x, temp(2), omega, alpha
	Integer i
	data omega /3.0d0/
	data alpha /0.5d0/
	If (i .EQ. 1) deriv=temp(2)
	If (i .EQ. 2) deriv= - omega**2 * temp(1) - alpha * temp(2)
Return
```
### rk4 

```fortran
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
```
### Result of x vs t
![](https://s31.postimg.org/g71fbtm1n/Screenshot_2016_07_18_21_30_26.png)

## Appendix for code
__P1__

```fortran
 	Program P1
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
        print *, sqrt(4.0D0**2.D0+2.d0*400*max/70)
c open file to output
        Open(6, File='p1.dat')
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
	Function deriv(temp, i)
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
	data alpha /0.5d0/
c
	If (i .EQ. 1) deriv=temp(2)
	If (i .EQ. 2) deriv=omega / temp(2)
	Return
	End

```

__P2__

```fortran
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
```

__P3__

```fortran
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
```
