set xlabel "Lattice Position" offset 0,0.5,0
set ylabel "Density" offset 2.75,0,0

plot "density-profile-A_1000000-iterations_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with linespoints pointtype 12 pointsize 1.5 lc rgb 'black' title "Component A Density Profile", \
 "density-profile-B_1000000-iterations_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with linespoints pointtype 6 pointsize 1.5 lc rgb 'red' title "Component B Density Profile"
