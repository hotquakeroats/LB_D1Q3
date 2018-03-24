#set terminal epslatex size 4.5,4.5 color solid
#set output "LB_nA00.800000_nB00.800000.eps"

set multiplot layout 3,1

set xlabel "Lattice Position" offset 0,0.5,0
set ylabel "Density" offset 2.75,0,0
set key outside center above
plot "density-profile-A_1000000-iterations_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with linespoints pointtype 12 pointsize 1.5 lc rgb 'black' title "Component A", \
 "density-profile-B_1000000-iterations_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with linespoints pointtype 6 pointsize 1.5 lc rgb 'red' title "Component B"

set xlabel "Lattice Position" offset 0,0.5,0
set ylabel "Pressure" offset 2.75,0,0
unset key
plot "pressure-profile_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with linespoints pointtype 12 pointsize 1.5 lc rgb 'black' title "Pressure", \
 "pressure-profile-theoretical_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'red' title "Theoretical Pressure"

set xlabel "Lattice Position" offset 0,0.5,0
set ylabel "Chemical Potential" offset 2.75,0,0
set format y '%e'
unset key
plot "chemical-potential-profile-A_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with linespoints pointtype 12 pointsize 1.5 lc rgb 'black' title "Component A Chemical Potential", \
 "chemical-potential-profile-B_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with linespoints pointtype 6 pointsize 1.5 lc rgb 'red' title "Component B Chemical Potential", \
 "chemical-potential-profile-theoretical-A_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' title "Theoretical Chemical Potential (A)", \
 "chemical-potential-profile-theoretical-B_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'red' title "Theoretical Chemical Potential (B)"

unset multiplot
