program P5
  implicit none
  real :: x1,x2,y1,y2,z1,z2,di
  real, external :: f

  print *,'please input the first point'
  read *,x1,y1,z1
  print *,'please input the second point'
  read *,x2,y2,z2

  call dist(x1,y1,z1,x2,y2,z2,di)
  print *, 'the distance by subroutine is ',di
  di = f(x1,y1,z1,x2,y2,z2)
  print *, 'the distance by function is ',di
end program P5

subroutine dist(x,y,z,a,b,c,d)
  implicit none
  real, intent(in) :: x,y,z,a,b,c
  real, intent(out) :: d
  d=sqrt((x-a)**2+(y-b)**2+(z-c)**2)
end subroutine dist

real function f(x,y,z,a,b,c)
  implicit none
  real, intent(in) :: x,y,z,a,b,c

  f=sqrt((x-a)**2+(y-b)**2+(z-c)**2)
end function f
