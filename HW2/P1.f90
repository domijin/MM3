program P1
  !  integrate exp(-x) from [-1,+1] by trapezoid algorithm
  implicit none
  integer :: i
  integer, dimension(4) :: N
  real*8 :: f_n,f_a,x0=-1,xf=1
  real*8, external :: f
  N = (/10,50,100,200/)
  f_a = -exp(x0)+exp(xf)
  
  do i = 1,4
     call trap(f,x0,xf,N(i),f_n)
     print *,N(i),f_n-f_a,(f_n-f_a)/f_a
  enddo
  
end program P1

subroutine trap(f,x0,xf,N,sum)
  implicit none
  real*8, intent(in) :: x0,xf
  integer, intent(in) :: N
  real*8, external :: f
  real*8, intent(out) :: sum
  real*8 :: h
  integer :: i
  sum = 0.5 * f(x0)
  h = (xf-x0)*1.0/N
  do i=1,N
     sum = sum + f(x0+i*h)
  enddo
  sum = (sum + f(xf))*h
end subroutine trap

real*8 function f(x)
  implicit none
  real*8, intent(in) :: x
  f=exp(-x)
end function f
