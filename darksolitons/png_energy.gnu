set terminal png
set output "energy.png"
set multiplot


set size 0.55,0.525
set origin -0.1,0.0
set title "POTENCIAL QUIMICO"
set yrange[41.8:42.0]
set xrange[0.0:30.0]
set xlabel "t"
set ylabel ""
plot "energies.txt" u 1:3 w l t""

set size 0.55,0.525
set origin 0.4,0.0
set title "POTENCIAL/INTERNA"
set yrange[7.5:20.0]
set xrange[0.0:30.0]
set xlabel "t"
set ylabel ""
plot "energies.txt" u 1:5 w l t"energia potencial" , "energies.txt" u 1:6 w l t"energia interna"

set size 0.55,0.525
set origin 0.4,0.5
set title "ENERGIA CINETICA"
set yrange[0.3:0.35]
set xrange[0.0:30.0]
set xlabel "t"
set ylabel ""
plot "energies.txt" u 1:4 w l t""

set size 0.55,0.525
set origin -0.1,0.5
set title "ENERGIA MEDIA"
set yrange[25.3:25.4]
set xrange[0.0:30.0]
set xlabel "t"
set ylabel ""
plot "energies.txt" u 1:2 w l t""

unset multiplot