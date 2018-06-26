        program generate_rand
        implicit none
        integer i

        do i=1, 10
            print *, rand()
C           if(mod(i,5).eq.1) print *, rand()
        enddo
            print *,'11th=', rand()
        end

