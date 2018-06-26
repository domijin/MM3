#P2-a: x-y path for the walk
set term jpeg size 800,1200;
set out "2-a.jpeg";
set multiplot
set autoscale

set origin 0,0.5
set size 1,0.5
set title "positive random walk"
plot "walk1.dat" every ::::499 u 1:2 notitle w lp lw 3

set origin 0,0
set size 1,0.5
set title "random random walk"
plot "walk2.dat" every ::::499 u 1:2 notitle w lp lw 3

unset multiplot;
unset out;

#P2-b: sqrt(R) vs. sqrt(N)
set term jpeg size 800,1200;
set out "2-b.jpeg";
set multiplot
set autoscale
set xlabel "sqrt(R)"
set ylabel "sqrt(N)"

set origin 0,0.5
set size 1,0.5
set title "positive random walk"
plot "walk1.dat" every ::::499 u 3:4 notitle 

set origin 0,0
set size 1,0.5
set title "random random walk"
plot "walk2.dat" every ::::499 u 3:4 notitle 

unset multiplot;
unset out;

#P2-c: sqrt(R) vs. sqrt(N) over 100 trials
set term jpeg size 800,1200;
set out "2-c.jpeg";
set multiplot
set autoscale
set xlabel "sqrt(R)"
set ylabel "sqrt(N)"

set origin 0,0.5
set size 1,0.5
set title "positive random walk"
plot "< awk '{R[$5]=R[$5]+$3;N[$5]=N[$5]+$4; nr[$5]++} END {for (i in R) {print i, R[i]/nr[i], N[i]/nr[i]}}' walk1.dat | sort -n" u 2:3 notitle 

set origin 0,0
set size 1,0.5
set title "random random walk"

plot "< awk '{R[$5]=R[$5]+$3;N[$5]=N[$5]+$4; nr[$5]++} END {for (i in R) {print i, R[i]/nr[i], N[i]/nr[i]}}' walk2.dat | sort -n" u 2:3 notitle 

unset multiplot;
unset out;