	program matrix
	parameter (n=100)
	real a(n,n), work(4*n)
        real gamma(n,n)
	real w(n)
	
	lwork = 4*n
c initialise
	do i=1,n
c diagonal element
	  a(i,i) = float(i)
 	  do j=i+1,n
c diagonal element
 	    a(i,j) = 1.0*float(i)/float(i+j)
	  enddo
	enddo

	call ssyev('V', 'U', n, a, n, w, work, lwork, info)
	print*,'info= ', info
	do i=1,10
  	  print*,w(i)
	enddo

C check eigenvector's ortho-normal
        call sgemm('T', 'N', n, n, n, 1.0, A, n, A, n, 0.0, gamma, n)
        do i=1,5
	  write(6,200) (gamma(j,i), j=1,10)
	enddo
200     format(1x,10f10.4)
	stop
	end
