program P2
  !  integrate exp(-x) from [-1,+1] by Simpson algorithm
  implicit none
  integer :: i
  integer, dimension(4) :: N
  real*8 :: f_n,f_a,x0=-1.,xf=1.
  real*8, external :: f
  N = (/10,50,100,200/)
  f_a = -exp(x0)+exp(xf)

  do i = 1,4
     call simp(f,x0,xf,N(i),f_n)
     print *,N(i),f_n-f_a,(f_n-f_a)/f_a
  enddo
  
end program P2

subroutine simp(f,x0,xf,N,sum)
  implicit none
  real*8, intent(in) :: x0,xf
  integer, intent(in) :: N
  real*8, external :: f
  real*8, intent(out) :: sum
  real*8 :: h
  integer :: i
  sum = -f(x0) + f(xf)
  h = (xf-x0)/2./N
  do i=0,2*N,2
     sum = sum + 4.*f(x0+h+i*h) + 2.*f(x0+i*h)
  enddo
  sum = sum*h/3.
end subroutine simp

real*8 function f(x)
  implicit none
  real*8, intent(in) :: x
  f=exp(-x)
end function f
