set xrange[-128:128]
set yrange[-0:0.3]
set xlabel "x"
set ylabel "|psi|^2"
set title ("Bright soliton")

plot "evolution_real.dat" u 1:2 index i with lines title "|psi(x)|^2",\
     "initial.dat" using 1:2 with line title "potential",\
     "initial.dat" using 1:3 with line title "initial state"

i=1+i
pause 0.005
if (i<201) reread
