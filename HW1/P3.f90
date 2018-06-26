program P3
  implicit none

  integer :: n,i,j,temp
  integer,dimension(:), allocatable :: serie
  print *,'please specify number of integers to sort'
  read *,n
  allocate(serie(n))
  print *,'please input',n,'integer line by line'
  do i=1,n
    read *,serie(i)
  enddo

  do i=1,n-1
    do j=2,n
        if (serie(i)<serie(j)) then
           temp=serie(i)
	   serie(i)=serie(j)
	   serie(j)=temp
	endif
    enddo
  enddo

  do i=1,n
   print *,serie(i)
  enddo
  deallocate(serie)
end program P3
