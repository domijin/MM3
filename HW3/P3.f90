program P3
  implicit none
  real*8, external :: drand48
  real*8 :: I=0,x
  integer :: n,j
  n = 0
  ! do-while loop until reach desired accuracy
  do while (abs(I-155.0/6)>1E-4)
     x = 0
     n = n + 1
     ! do loop to scatter one point in 10-D space
     do j = 1,10
       x = x + drand48()
     enddo
     I = I * ( n - 1 ) + x**2 ! reverse the total value in 10-D space for n trails
     I = I / n ! average value of all the n possible points is the integration result
     print "(i6,2x,f6.3)", n, I
  enddo
end program P3