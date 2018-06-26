# E & M vs T

set term jpeg size 800,1200;
set out "Result.jpeg";
set multiplot
set autoscale

set origin 0,0.5
set size 1,0.5
set title "Magnetization aganist T"
plot "res.dat" u 5:3 notitle w lp lw 3

set origin 0,0
set size 1,0.5
set title "Energy aganist T"
plot "res.dat" u 5:1 notitle w lp lw 3

unset multiplot;
unset out;
