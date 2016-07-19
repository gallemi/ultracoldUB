set multiplot

set title "TRAYECTORIA SOLITON"
set yrange[0.0:0.1]
set xrange[-15.0:15.0]
set xlabel "z"
set ylabel "abs(psi)**2"
plot "WfDs.txt" u 1:2 index i with lines t "density"

i=1+i

pause 0.1
if (i<1200) reread
unset multiplot

