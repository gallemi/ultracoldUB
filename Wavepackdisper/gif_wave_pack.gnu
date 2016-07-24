do for [l=1:1000]{
k=sprintf('%08d.txt',l)
set multiplot layout 2,2


set title "DISPERSION PAQUETE DE ONDAS"
set yrange[0:1.2]
set xrange[-35.0:35.0]
set xlabel "z"
set ylabel "abs(psi)**2"
plot 'WfWd-'.k using 1:2 w l t""



set title "DISPERSION PAQUETE DE ONDAS"
set yrange[-1.0:1.0]
set xrange[-35.0:35.0]
set xlabel "z"
set ylabel "psi"
plot 'WfWd-'.k using 1:4 w l t"real wf"


set title "DISPERSION PAQUETE DE ONDAS"
set yrange[0.0:0.1]
set xrange[-35.0:35.0]
set xlabel "z"
set ylabel "abs(psi)**2"
plot 'WfWd-'.k using 1:6 w l t"K-Space"


set title "DISPERSION PAQUETE DE ONDAS"
set yrange[-1.0:1.0]
set xrange[-35.0:35.0]
set xlabel "z"
set ylabel "psi"
plot 'WfWd-'.k using 1:5 w l t"imag wf"; pause 0.2}

unset multiplot


