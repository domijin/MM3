      real*8 function drand48()
!-------------------------------------------------------------------------
!
      integer m, ia, ic, ntab
      real*8  rm
      parameter (ntab=97,m=714025,ia=1366,ic=150889,rm=1.0/m)
      integer ir(ntab), iff, idum, j, iy
      data iff /0/, idum/0/
      save iff, idum, iy, ir
!
!
      if(idum.lt.0.or.iff.eq.0) then
        iff=1
        idum=mod(ic-idum,m)
        do j=1,ntab
           idum=mod(ia*idum+ic,m)
           ir(j)=idum
        end do
        idum=mod(ia*idum+ic,m)
        iy=idum
      endif
      j=1+(ntab*iy)/m
      if(j.gt.ntab.or.j.lt.1) then
        print*, 'error! in drand48, j out of range',j
	  stop
      endif
      iy=ir(j)
      drand48=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end

