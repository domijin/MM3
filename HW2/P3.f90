program P3
  !  integrate exp(-x) from [-1,+1] by Gaussian algorithm
  implicit none
  integer :: i,j
  real*4, dimension(:), allocatable :: x,w
  integer, dimension(2) :: n
  real*4 :: f_n,f_a,x0=-1,xf=1
  real*4, external :: f

  n = (/4,8/)
  f_a = -exp(x0)+exp(xf)
  
  do i = 1,2
     allocate(x(n(i)),w(n(i)))
     call gauleg(x0,xf,x,w,n(i))

     f_n = 0
     do j = 1,n(i)
        f_n = f_n + w(j)*f(x(j))
     enddo
     
     print *,n(i),f_n-f_a,(f_n-f_a)/f_a
     deallocate(x,w)
  enddo
  
end program P3

real*4 function f(x)
  implicit none
  real*4, intent(in) :: x
  f=exp(-x)
end function f
