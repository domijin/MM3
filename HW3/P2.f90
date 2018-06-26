program P2
  implicit none
  integer :: i
  open(unit=20,file='walk1.dat')
  open(unit=21,file='walk2.dat')
  do i=1,100
    call walk1
    call walk2
  enddo
  close(20)
  close(21)
  call system("gnuplot plot.gnu")
  print *,"The distance random walk can reach should be propotional to the sqrt(steps). Plot 2-c.jpeg agrees with my conclusion."
  ! plot.gnu will generate 2-[a,b,c].jpeg, 2-[a,b].jpeg will only use first 500 lines in data files
end program P2

subroutine walk1
  ! this version will only walk to x+ and y+ directions
  implicit none
  real*8, external :: drand48
  integer :: n, x, y
  x=0
  y=0
  do n=1,500
    ! map rand to [-1, 1], positive value for heading to y direction; negative value for heading to x direction.
     if (2*drand48()-1.le.0) then
        x=x+1
     else
        y=y+1
     endif
     write (20,*),x,y,sqrt(x*x*1.0+y*y*1.0),sqrt(n*1.0),n
     ! output x,y,sqrt(R),sqrt(N),n to plot.dat for plotting
  enddo
end subroutine walk1

subroutine walk2
  ! this version will walk to x+,x-,y+,y- directions
  implicit none
  real*8, external :: drand48
  integer :: n, x, y
  x=0
  y=0
  do n=1,500
    ! map rand to [-1, 1], positive value for heading to y direction; negative value for heading to x direction.
     if (2*drand48()-1.le.0) then
        ! map rand to [-1, 1], positive value for heading to x+ direction; negative value for heading to x- direction.
        if (2*drand48()-1.le.0) then
          x=x-1
        else
          x=x+1
        endif
     else
        ! map rand to [-1, 1], positive value for heading to y+ direction; negative value for heading to y- direction.
        if (2*drand48()-1.le.0) then
          y=y-1
        else
          y=y+1
        endif
     endif
     write (21,*),x,y,sqrt(x*x*1.0+y*y*1.0),sqrt(n*1.0),n
     ! output x,y,sqrt(R),sqrt(N),n to plot.dat for plotting
  enddo
end subroutine walk2