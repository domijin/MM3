program p1

  implicit none
  real :: F,C

  print *,'Please input a emperature in F'
  read *,F
  open(unit=20,file='output.dat')

  C=(F-32)*5./9
  print *,C
  write (20,*) F,' F',achar(10),C,' C'
  close(20)

end program
