set multiplot

set size 0.45,0.50
set origin 0.0,0.0
set title "POSICION DEL MINIMO DE DENSIDAD"
set yrange[-4.0:4.0]
set xrange[0.0:35.0]
set xlabel "t"
set ylabel "position"
plot "min.txt" u 1:2 w l t"" , "min.txt" every::i::i u 1:2 points 1 pointtype 8 t""

set size 0.45,0.50
set origin 0.0,0.5
set title "VALOR DEL MINIMO DE DENSIDAD"
set yrange[-0.1:0.1]
set xrange[0.0:35.0]
set xlabel "t"
set ylabel "value"
plot "min.txt" u 1:3 w l t"" , "min.txt" every::i::i u 1:3 points 1 pointtype 8 t""

set size 0.45,0.50
set origin 0.5,0.5
set title "ENERGIA CINETICA"
set yrange[0.3:0.35]
set xrange[0.0:35.0]
set xlabel "t"
set ylabel "energia"
plot "energies.txt" u 1:4 w l t"" , "energies.txt" every::i::i u 1:4 points 1 pointtype 8 t""

set size 0.45,0.50
set origin 0.5,0.0
set title "ENERGIA POTENCIAL"
set yrange[7.5:10.0]
set xrange[0.0:35.0]
set xlabel "t"
set ylabel "energia"
plot "energies.txt" u 1:5 w l t"" , "energies.txt" every::i::i u 1:5 points 1 pointtype 8 t""

i=1+i

pause 0.1
if (i<340) reread
unset multiplot