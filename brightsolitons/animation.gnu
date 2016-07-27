do for [l=0:20000:100]{
k=sprintf('%08d.dat',l)

a=l/1000
set title sprintf("Bright soliton (t=%2.1f)",a)
#set title ("Bright soliton")
set xrange[-128:128]
set yrange[-0:0.3]
set xlabel "x"
set ylabel "|psi|^2"
plot './bs_evolution/WfBs-'.k using 1:2 w l t "|psi(x)^2|",\
	 'initial.dat' using 1:2 w l t "potential",\
	 'initial.dat' using 1:3 w l t "initial state"

}
