program ising_model_2d
	implicit none
	real*8, external :: drand48
	integer :: i,j,x,y,nmc=10000,n,transient
	real*8 :: T,e,M,dE, Etot=0,Mtot=0, E_avg=0,M_avg=0,norm,temp_e
	integer, dimension(2) :: lat(20,20)

	norm =1.d0/4.0D6
	! initialize the lattice
	! python: lat.append(random.choice([1,-1],size=L*L))
	do x=1,20
		do y=1,20
			if (drand48().ge.0.5) then 
				lat(x,y)=1
			else
				lat(x,y)=-1
			end if
		enddo
	enddo

	do T=5.0,0.0,-0.1
		do transient=1,999 
			do n=1,400
				! choose random pos
				x=floor(drand48()*20)+1
				y=floor(drand48()*20)+1
				! check energy
				! merge function is fortran version of lambda function. merge(resA,resB,cond), if cond true, resA, else resB
				e=-1*lat(x,y)*(lat(merge(20,x-1,x<2),y)+lat(merge(1,x+1,x>19),y)+lat(x,merge(1,y+1,y>19))+lat(x,merge(20,y-1,y<2)))
				dE=-2*e
				if (dE.lt.0 .or. drand48().le.exp(-dE/T)) then
					! flip
					lat(x,y)=-lat(x,y)
				end if
			enddo
		enddo
		! ceil energy calculation in python:
		!	E=-1*sum(multiply(lat,roll(lat,1,axis=0)+roll(lat,-1,axis=0)+roll(lat,1,axis=1)+roll(lat,-1,axis=1)))

		M=0
		E=0

		! total magnetization
		do x=1,20
			do y=1,20
				M=M+lat(x,y)
				E=E-1*lat(x,y)*(lat(merge(20,x-1,x<2),y)+lat(merge(1,x+1,x>19),y)+lat(x,merge(1,y+1,y>19))+lat(x,merge(20,y-1,y<2)))
			enddo
		enddo
		Etot=0
		Mtot=0

		do i=1,nmc ! Monte Carlo loop
			do j=1,400 ! Metropolis loop
				! choose random pos
				x=floor(drand48()*20)+1
				y=floor(drand48()*20)+1
				! check energy
				temp_e=-1*lat(x,y)*(lat(merge(20,x-1,x<2),y)+lat(merge(1,x+1,x>19),y)+lat(x,merge(1,y+1,y>19))+lat(x,merge(20,y-1,y<2)))
				dE=-2*temp_e
				if (dE.lt.0 .or. drand48().le.exp(-dE/T)) then
					! flip
					lat(x,y)=-lat(x,y)
					! update E, M
					E=E+2*dE
					M=M+2*lat(x,y)
				endif
			enddo
			Etot=Etot+E/2.0
			Mtot=Mtot+M
		enddo
		E_avg=Etot*norm
		M_avg=Mtot*norm;
		print *, E,E_avg,M,M_avg,T
	enddo
end program

