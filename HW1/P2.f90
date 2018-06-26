program P2
  implicit none
  real :: x1,y1,z1,x2,y2,z2
  print *,'input the first point'
  read *,x1,y1,z1
  print *,'input the second point'
  read *,x2,y2,z2
  print *,'distance is',sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

end program P2
