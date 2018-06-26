        program generate_rand
        implicit none
        integer i
        real rand
	real*8  drand48

        do i=1, 20
            print *, rand(), drand48()
        enddo
        end
