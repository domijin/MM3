      	Program integrate
      	Implicit none
c declarations 
      	Real*8 trapez, simpson, quad, r1, r2, r3
      	Real*8 theo, vmin, vmax
      	Integer i,j, id(8)
	data id/2, 10, 20, 40, 80, 160, 320, 640/
c
c theoretical result, integration range
      	theo = 0.632120558829
c	theo = dexp(1.0D0)
c	theo = theo - 1.d0/theo
      	vmin=0.0
      	vmax=1.0
c     	Open(6, File='integ.dat', Status='Unknown')
c calculate integral using both methods for steps = 3..501
        Do 50 j=1, 8
	  i = id(j) + 1
          r1=trapez(i, vmin, vmax)
          r1=abs(r1-theo)
          r2=simpson(i,vmin, vmax)
          r2=abs(r2-theo)
          r3=quad(i,vmin, vmax)
          r3=abs(r3-theo)
C         write(6,100) i, log10(r1), log10(r2), log10(r3)      
          write(6,100) i, r1, r2, r3      
 50   Continue
 	Close(6)
      	Stop 'data saved in integ.dat'
C100     format(1x,i7,1x,3f12.5)
100     format(1x,i7,1x,3e13.4)
      	End
c
c the function we want to integrate
      	Function f(x)
      	Implicit none
      	Real*8 f, x
          f=dexp(-x)
       	Return
      	End
c
c trapezoid rule
      	Function trapez(i, min, max)
      	Implicit none
      	Integer i, n	
      	Real*8 f, interval, min, max, trapez, x
      	trapez=0 
      	interval = ((max-min) / (i-1))
c sum the midpoints
      	Do 21 n=2, (i-1)   
          x = interval * (n-1) + min
          trapez = trapez + f(x)*interval
 21   	Continue 
c add the endpoints  
      	trapez = trapez+0.5*(f(min)+f(max))*interval
        Return
      	End
c
c Simpson's rule
      	Function simpson(i, min, max)
      	Implicit none
      	Integer i, n
      	Real*8 f, interval, min, max, simpson, x
      	simpson=0	
      	interval = ((max-min) / (i-1))
c loop for odd points
      	Do 31 n=2, (i-1), 2        
          x = interval * (n-1) + min
          simpson = simpson + 4*f(x)
 31   	Continue 
c loop for even points
      	Do 32 n=3, (i-1), 2          
          x = interval * (n-1) + min
          simpson = simpson + 2*f(x)
 32   	Continue  
c add the endpoints  
      	simpson = simpson+f(min)+f(max)
      	simpson=simpson*interval/3
        Return
      	End
c
c Gauss' rule
      	Function quad(i, min, max)
      	Implicit none
      	Real*8 w(1000), x(1000)
      	Real*8 f, min, max, quad
      	Integer i, job, n
      	quad=0
      	job=0
      	call gauss(i, job, min, max, x, w)
      	Do 41 n=1, i
          quad=quad+f(x(n))*w(n)
 41   	Continue
      	Return 
      	End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  gauss.f: Points and weights for Gaussian quadrature                 c
c								       c
c  taken from: "Projects in Computational Physics" by Landau and Paez  c
c	       copyrighted by John Wiley and Sons, New York            c
c                                                                      c
c  written by: Oregon State University Nuclear Theory Group            c
c	       Guangliang He & Rubin H. Landau                         c
c  supported by: US National Science Foundation, Northwest Alliance    c
c                for Computational Science and Engineering (NACSE),    c
c                US Department of Energy 	                       c
c								       c
c  comment: error message occurs if subroutine called without a main   c
c  comment: this file has to reside in the same directory as integ.c   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     rescale rescales the gauss-legendre grid points and weights
c
c     npts     number of points
c     job = 0  rescalling uniformly between (a,b)
c           1  for integral (0,b) with 50% points inside (0, ab/(a+b))
c           2  for integral (a,inf) with 50% inside (a,b+2a)
c     x, w     output grid points and weights.
c
      subroutine gauss(npts,job,a,b,x,w) 
      integer npts,job,m,i,j 
      real*8 x(npts),w(npts),a,b,xi
      real*8 t,t1,pp,p1,p2,p3,aj
      real*8 eps,pi,zero,two,one,half,quarter
      parameter (pi = 3.14159265358979323846264338328, eps = 3.0E-14)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      parameter (half=0.5d0,quarter=0.25d0)
c
c *** FIRST EXECTUABLE *************************************************
c

      m=(npts+1)/2
      do 1020 i=1,m
         t=cos(pi*(i-quarter)/(npts+half))
 1000    continue
         p1=one
         p2=zero
         aj=zero
         do 1010 j=1,npts
            p3=p2
            p2=p1
            aj=aj+one
            p1=((two*aj-one)*t*p2-(aj-one)*p3)/aj
 1010    continue
         pp=npts*(t*p1-p2)/(t*t-one)
         t1=t
         t=t1-p1/pp
c
         if(abs(t-t1).gt.eps) goto 1000
c
         x(i)=-t
         x(npts+1-i)=t
         w(i)=two/((one-t*t)*pp*pp)
         w(npts+1-i)=w(i)
 1020 continue
c
c rescale the grid points 
      if (job.eq.0) then
c     scale to (a,b) uniformly
         do 1030 i=1,npts
            x(i)=x(i)*(b-a)/two+(b+a)/two
            w(i)=w(i)*(b-a)/two
 1030    continue
      elseif (job.eq.1) then
c scale to (0,b) with 50% points inside (0,ab/(a+b))
         do 1040 i=1,npts
            xi=x(i)
            x(i)=a*b*(one+xi)/(b+a-(b-a)*xi)
            w(i)=w(i)*two*a*b*b/((b+a-(b-a)*xi)*(b+a-(b-a)*xi))
 1040    continue
      elseif (job.eq.2) then
c scale to (a,inf) with 50% points inside (a,b+2a)
         do 1050 i=1,npts
            xi=x(i)
            x(i)=(b*xi+b+a+a)/(one-xi)
            w(i)=w(i)*two*(a+b)/((one-xi)*(one-xi))
 1050    continue
      else
         pause 'Wrong value of job'
      endif
c
      return
      end
