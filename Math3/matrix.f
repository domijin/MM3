	program matrix
	parameter (n=700)
	real a(n,n), b(n,n), c(n,n)
	
c initialise
	do i=1,n
	  do j=1,n
	    a(i,j) = 1.0/float(i+j)
	    b(i,j) = 1.0/float(i) + 1.0/float(j)
	  enddo
	enddo

	do i=1,n
	  do j=1,n
	    c(i,j) = 0.0
	    do k=1,n
	      c(i,j) = c(i,j) + a(i,k)*b(k,j)
	    enddo
	  enddo
	enddo	
	print*,c(10,10)
	stop
	end
