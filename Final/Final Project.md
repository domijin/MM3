# Final Project
Path to source code: `/home/d/dx/dxj4360/Final`

File list:

`final.f90`: scan method to find the lowest 4 eigenvalues  
`final2.f90`: bisection method with advancing E to find the lowest 4 eigenvalues  
`manual.f90`: use bisection method to find eigenvalue between user input e1,e2  
`psive.f90`: loop over E to find the tendency of $\psi(x_0)$ for eigenvalues

## The Problem
The one-dimensional double-well problem:

$[-\frac{1}{2}\frac{d^2}{dx^2}+V(x)]\psi(x)=E \psi(x)$

with potential $V(x)=ax^4-bx^2$, where $a=2,b=8$

[](https://s10.postimg.org/64qwjilop/Screenshot_2016_08_08_15_32_47.png)

### The ODE
It can be reformulated to:

$\psi''=2(V(x)-E)\psi=(4x^4-16x^2-2E)\psi$

which doesn't include $\psi'$. Thus, 

$\frac{d\psi}{dx}=\psi'$

$\frac{d\psi'}{dx}=(4x^4-16x^2-2E)\psi$

### The boundary condition & initial parameters

It is impossible to integrate the wave function over $[-\infty,\infty]$ and archive unity. But the potential grows very fast outside $|x|>x_0$, where $x_0$ is set to 2.5 in my calculation. Thus, the problem is approximately solved as a double-well potential within infinite potential box.

* lower bound energy level: $E=\min(V(x))=-8$
* infinity assumption: $x_0=2.5$ such that $\lim_{|x|\ge x_0}\psi(x)\rightarrow\epsilon$; $\psi'(x_0)\sim \epsilon$
* resolution: $dE=1e-4, dx=1e-3$
* $\psi(-x_0)=\psi'(-x_0)=\epsilon$
* converge condition: $\psi(x_0)\le \epsilon$

## Pseudocode

```Fortran
    initialize parameters
    do loop over E
        runge kutta integration from -x_0 to x_0
        if psi(x_0) -> 0 then
            record all psi(x)
            print E
        E increase by some steps         
```

## Result

* (a) energy eigenvalue below 0 when $x_0=2.5$:
    * -5.3229
    * -5.3035 
    * -1.0206
    * -0.3616
* (b) lowest two states with $V(x)$
 ![](https://s9.postimg.org/q7qjmulsv/Screenshot_2016_08_09_21_39_21.png)

* (c) if $V(x)$ is approximated as SHO at $x_0=\pm\sqrt2$, the wave function will be symmetric at $x=x_0$, the energy eigenvalue will change depend on the shape of SHO potential. If the minimum potential is the same, the wider potential, the higher the energy eigenstate.
![](https://s9.postimg.org/yf7zm7qm7/Screenshot_2016_08_09_14_01_03.png)


## Discussion

Manually guessing the energy eigenvalue is very inefficient, so I used do-loop to scan from $E=[-8,0]$ with some small steps. I found the eigenstates are highly dependent on the shooting interval. If $x_0=\pm3$ is assumed as the cut-off boundary, the ground state energy will be extremely hard to find, even with a more efficient bisection method and adapting dx resolution. So I checked the $\psi(x_0)$ vs $E$. The eigenvalues are extremely narrow and easy to miss. If $dx$ is too small, $\psi$ will accumulate a lot at the beginning and hard to cross 0. 