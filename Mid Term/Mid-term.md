## Mid-term 
by Dongming Jin

###I multiple choices
1. c 
2. b
3. c
4. d
5. b
6. b
7. d
8. c
9. c
10. b

### II

* a)

$\int_a^bf(x)dx \\ 
\simeq (\frac{1}{2}(f_1+f_2)+\frac{1}{2}(f_2+f_3)+...+\frac{1}{2}(f_N+f_{N+1})) * \frac{b-a}{N} \\
=\frac{h}{2}\Sigma_{k=1}^{N+1}(f(x_{k+1})+f(x_k))\\
=h(\frac{1}{2}f_1+f_2+...+f_N+\frac{1}{2}f_{N+1})$  


where $h=(b-a)/N$

* b) 

$\int_{x_i}^{x_{i+1}}f(x)dx \simeq \int_{x_i}^{x_{i+1}}[f_i+(x-x_i)f'_i+\frac{1}{2}(x-x_i)^2f_i'']dx$


_I get stuck for the expansion and forget the trick. It should be derived by Tyler expansion minus the trapezoid equation. The result I remember is about $\frac{1}{12}h^3f''$_

###III

* a)

$\frac{dx}{dt}=v$

$\frac{dv}{dt} = \frac{d^2x}{xt^2}= -\omega^2 x -\alpha v$

* b)

$x_{i+1}-x_i=hv_i$ ->> $x_{i+1} = x_i+hv_i$

$v_{i+1}-v_i= h(-\omega^2 x_i - \alpha v_i)$ ->> $v_{i+1} = v_i -h\omega^2 x_i - h\alpha v_i$

### IV

* a) 

$\frac{\text{area of shadow}}{\text{area of square}}=\frac{1}{4}\pi r^2= \frac{\text{points in shadow}}{\text{total points}}$

$\iint dxdy=\frac{\text{points in shadow}}{\text{total points}}$

__procedure__

```fortran
do i = 1,max
    x=drand48()
    y=drand48()
    if (x*x+y*y.le.1 ) pts = pts + 1
enddo
print *, 'integral is ', pts/max
```

* b)

```fortran
do i = 1,max
    x=drand48()
    value = value + f(x)
enddo
print *, 'integral is ', value/max
```

###V

The central method is more accurate.

* the forward method:
$f(x+h)\simeq f(x)+hf'(x)+\frac{1}{2}h^2f''(x)+ O(h^3f''')+...$ (1)
thus $f'_{fd}(x)\simeq\frac{f(x+h)-f(x)}{h}-\frac{1}{2}hf''(x)\sim O(hf'')$

* the central method
$f(x-h)\simeq f(x)-hf'(x)+\frac{1}{2}h^2f''(x)-O(h^3f''')+...$ (2)   
(1)-(2): $f(x+h)-f(x-h)=2hf'(x)+2O(h^3f''')$  
$f'_{cd}(x)\simeq \frac{f(x+h)-f(x-h)}{2h}-O(h^2f''')\sim O(h^2f''')$

###VI

No. Computer has finite accuracy and limited computing power. For double precision float, it use 53bit for decimal points, which means the minimum it can represent is $2^{-53} \simeq10^{-16}$. Below that, the truncation error will domain and accuracy won't increase anymore. 