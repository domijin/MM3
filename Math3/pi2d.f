	Program pai
	Implicit none
	integer max
c declarations
	Real volume, x, y, area
	Integer i, pi3, pi2

	print*,'input the number for try'
	read(5,*) max
	pi2=0
c execute
	Do 10 i=1, max
C generate (x,y) within [-1,1]:
	   x = rand()*2-1
	   y = rand()*2-1	
	   If ((x*x + y*y ) .LE. 1) pi2 = pi2 + 1
	   area = 4.0 * pi2/Real(i)

 	   if(mod(i,1000).eq.1) Write(6,100) i, area
 10 	Continue

        Write(6,100) max, area
100	format(1x,'try= ',i8,2x,'pai=',f8.4)
	stop
	end
 	
