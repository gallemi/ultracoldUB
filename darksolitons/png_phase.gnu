set terminal png
set output "dif_phase.png"

set title "DIFERENCIA DE FASE DEBIDA A SOLITON"
set yrange[0.0:5.0]
set xrange[0.0:30.0]
set xlabel "t"
set ylabel "fase"
plot "phase.txt" u 1:2 t""