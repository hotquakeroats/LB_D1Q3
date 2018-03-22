set loadpath system("echo $LB_HOME/Maxwell Construction/single component/phase-data-4/")
set lmargin at screen 0.05
set rmargin at screen 0.95
set bmargin at screen 0.05
set tmargin at screen 0.95
set xlabel "Density"
set ylabel "Temperature Ratio (T / Tc)"
set grid

#set terminal epslatex size 4.5,4.5 color solid
#set output "singlecomponent-phase-diagram.eps"

plot "phaseDiagram_rhoVsTemp_theoretical.dat" with lines, \
 "phase-data-gradP-CC-periodic.dat", \
 "phase-data-gradP-AB-periodic.dat", \
 "phase-data-gradMuFull-CC-periodic.dat", \
 "phase-data-gradMuFull-AB-periodic.dat", \
 "phase-data-gradP-CC-walls.dat", \
 "phase-data-gradP-AB-walls.dat", \
 "phase-data-gradMuFull-CC-walls.dat", \
 "phase-data-gradMuFull-AB-walls.dat"

#set output
set terminal wxt persist
unset lmargin
unset rmargin
unset bmargin
unset tmargin
replot

