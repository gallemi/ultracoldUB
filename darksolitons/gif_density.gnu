do for [l=1:350]{
k=sprintf('%08d.txt',l)
set multiplot layout 2,2


set title "TRAYECTORIA SOLITON"
set yrange[0:0.1]
set xrange[-15.0:15.0]
set xlabel "z"
set ylabel "abs(psi)**2"
plot 'WfDs-'.k using 1:2 w l t""



set title "TRAYECTORIA SOLITON"
set yrange[-0.5:0.5]
set xrange[-15.0:15.0]
set xlabel "z"
set ylabel "psi"
plot 'WfDs-'.k using 1:4 w l t"real wf"


set title "TRAYECTORIA SOLITON"
set yrange[0.0:110.0]
set xrange[-15.0:15.0]
set xlabel "z"
set ylabel "psi"
plot 'WfDs-'.k using 1:6 w l t"potential"


set title "TRAYECTORIA SOLITON"
set yrange[-0.5:0.5]
set xrange[-15.0:15.0]
set xlabel "z"
set ylabel "psi"
plot 'WfDs-'.k using 1:5 w l t"imag wf"; pause 0.2}

unset multiplot


