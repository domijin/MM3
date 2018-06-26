	Program pai
	Implicit none
	integer max
c declarations
	Real volume, x, y, z, area
	Integer i, pi3, pi2
	print*,'input the number for try'
	read(5,*) max
	pi3=0
c use max as random seed
	call srand(max)
c execute
	Do 10 i=1, max
c generate (x,y) within [-1,1]:
	   x = rand()*2-1
	   y = rand()*2-1
	   z = rand()*2-1
	   If ((x*x + y*y + z*z) .LE. 1) pi3 = pi3 + 1
	   volume = 6.0 * pi3/Real(i) 
c volume inside sphere ~ pts/i = 4/3*pi*r**3 --> pi= pts*8
 	   if(mod(i,int(sqrt(i*1.0))).eq.1) Write(6,100) i, volume
 10 	Continue

        Write(6,100) max, volume
100	format(1x,'try= ',i10,2x,'pai=',f8.4)
c do-while loop to get enough significant numbers
	pi3=0
	i=0
	do while(abs(volume-3.1415927).GT.1E-4)
	   x = rand()*2-1
	   y = rand()*2-1
	   z = rand()*2-1
	   i = i+1
	   if ((x*x + y*y + z*z).LE.1) pi3 = pi3 +1
           volume = 6.0 * pi3/Real(i)
 	enddo
	write (6,200), i,volume
200	format(i8," trials are needed, the final value of Pi is ",f5.3)
	stop
	end
 	
