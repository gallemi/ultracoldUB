v = 6.0
gn = -0.4
h = 0.5
w = 0.5 

set xrange[-128:128]
set yrange[-0:0.3]
set xlabel "x"
set ylabel "|psi|^2"
set title sprintf("Bright soliton (v=%.1f, gn=%.1f, h/E=%.1f, w/xi=%.1f)", v, gn, h, w)

plot "evolution_real.dat" u 1:2 index i with lines title "|psi(x)|^2",\
     "initial.dat" using 1:2 with line title "potential",\
     "initial.dat" using 1:3 with line title "initial state"

i=1+i
pause 0.005
if (i<201) reread
