	parameter (NGX=256)
	real*8 x, w, a, h
        complex*16 f(NGX)	
	integer nn
	integer i, j

	nn = ngx
c stepsize
        h = 1.d0/ngx	
	x=0.d0
	do i=1, NGX
	  f(i) = dcmplx(dcos(6.28*x)*dcos(62.8*x), 0.d0)
	  x = x + h
	enddo
        do i=1,40,4
          print*,i,f(i)
	enddo
C FFT
	call four1( f, nn,-1)
C inv-FFT
	call four1( f, nn, 1)
        do i=1,40,4
c	  a = f(i)*conjg(f(i))
	  print*,i,f(i)/nn
	enddo
	stop
	end
	  
