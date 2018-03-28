set loadpath system("echo $LB_HOME/Maxwell Construction/single component/phase-data-4/")
set lmargin at screen 0.05
set rmargin at screen 0.95
set bmargin at screen 0.05
set tmargin at screen 0.95
set xlabel "Phase Density"
set ylabel "Temperature Ratio (T / Tc)"
set xrange [-0.1:3.1]
set yrange [0.0:1.1]
unset grid

#set terminal epslatex size 4.5,4.5 color solid
#set output "singlecomponent-phase-diagram.eps"

plot "phaseDiagram_rhoVsTemp_theoretical.dat" with lines linecolor rgb 'black' title "Binodal", \
 "phase-data.dat" with points pointtype 7 linecolor rgb 'blue' title "LB densities"

#set output
set terminal wxt persist
unset lmargin
unset rmargin
unset bmargin
unset tmargin
replot

