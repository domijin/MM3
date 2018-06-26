	program matrix_lib
	parameter (n=700)
	real a(n,n), b(n,n), c(n,n)
	
c initialise
	do i=1,n
	  do j=1,n
	    a(i,j) = 1.0/float(i+j)
	    b(i,j) = 1.0/float(i) + 1.0/float(j)
	  enddo
	enddo

	call sgemm('n', 'n', n, n, n, 1.0, a, n, b, n, 0.0, c, n)
	print*,c(10,10)
	stop
	end
