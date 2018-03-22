set terminal wxt font "Verdana,14"
set xlabel "Density"
set ylabel "Temperature"
set xrange [-0.1:3.1]
set yrange [0:1.1]
unset xtics
unset ytics

plot "phaseDiagram_rhoVsTemp_theoretical.dat" using 1:2 with lines lc rgb 'black' title "Binodal", \
 "phaseDiagram_rho-rhoVsTemp_tie-lines.dat" using 1:2 with linespoints pointtype 7 lc rgb 'red' title "Tie Lines"
