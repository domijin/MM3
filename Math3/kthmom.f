        program generate_rand
        implicit none
        integer i,k,n
	real sum
        real*8 drand48

	sum =0.0
	print*,'input k and N'
	read*,k,n
        do i=1, n
            sum = sum + drand48()**k
        enddo
	sum=sum/n
            print *,'the k-th moment =', sum
        end

