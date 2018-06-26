program P4
  implicit none
  real*4 :: pi,one
  real*8 :: pi8,one8

  one=1.
  one8=1.
  pi=atan(one)*4
  pi8=atan(one8)*4

  print "('pi in single ',f20.18,' and pi in double ',f20.18)", pi,pi8

end program P4

