# Project 1

Path to the source code: `/home/d/dx/dxj4360/Project1`
## Introduction
__Statistical Background__

* `Ferromagnetic` materials will tend to have ordered magnetic dipole momentum when exposed to external magnetic field and keep `ordered phase` afterwards. 
* `Paramagnetic` materials can be weakly induced by external magnetic field but will revert to `disordered phase` when the external field is removed. 
* `Phase transition` is when the microstates of materials change between ordered phase and disordered phase.
* A `microstate` is a specific microscopic configuration of a thermodynamic system. 
* Due to thermal fluctuations, each mircostate has a certain `probability` of occurrence, which is the possibility of the microscopic configuration: $P(x)=\frac{1}{Z}e^{-\beta H(x)}$. The total probability of all configurations is 1.
* `Statistical average`, in my understanding, is the ensemble of all possible states, which equals the mean value of all microstates of a system. 


__Ising Model__

* The `Ising Model` is a mathematical model in statistical mechanics. It uses discrete variables to represent magnetic dipole moments of atomic `spins`, whose `possible value` is $\pm1$. 
* The `energy` of a configuration $\sigma$ is given by the Hamiltonian function:
${\displaystyle H(\sigma )=-J\sum _{\langle ij\rangle }\sigma _{i}\sigma _{j}-h\sum _{j}\sigma _{j}.}$
* For a system of $N$ spins $(N=L\times L)$, there are $2^N$ microstates.

__Monte Carlo method__

* For a 2-dimensional square lattice at the no external magnetic field case ($H=0$), we have
    * $L = 20$: the total number of sites on the lattice,
    * $\sigma_j \in \{−1, +1\}$: an individual spin site on the lattice, $j = 1, ..., L$,
    * $S \in \{−1, +1\}^L$: state of the system.
* Since there is no external field, the Hamiltonian function is thus, 
${\displaystyle H(\sigma )=-J\sum _{\langle ij\rangle }\sigma _{i}\sigma _{j}.}$

* The probability function has an actual statistical weight. In a discrete case that the phase space a computer algorithm will generate is finite, we need to use `importance sampling` for better estimating the properties of a particular distribution.

* The total energy $E_{flip}$ can be calculated from the Hamiltonian given earlier: $\langle E \rangle=\frac{1}{2} \langle -J\sum _{\langle ij\rangle }\sigma _{i}\sigma _{j} \rangle$



## Pseudocode
![Flow Chart](https://s32.postimg.org/4le7ln2px/Flow_Chart.png)

## Fortran code
```fortran
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
		do transient=1,1000 
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
```

## Result

![Result](https://s31.postimg.org/3oy5hhwyh/Result.jpg)

The Curry temperature is around $2.2 J/k_B$

Above the critical temperature, the spontaneous magnetization vanishes as the thermo effect surpass the ferromagnetic state. 

From the first image, we can see the total Magnetization drops as the temperature increases. The theromo fluctuation of temperature increase starts to destroy the configuration. After reaching the critical temperature, the configuration is totally randomized, thus there is no preferred magnetic direction. Correspondingly, we can see the total energy from the second image increases as the temperature goes up.