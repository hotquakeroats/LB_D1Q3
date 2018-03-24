set xlabel "Lattice Position" offset 0,0.5,0
set ylabel "Pressure" offset 2.75,0,0

plot "pressure-profile_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with linespoints pointtype 12 pointsize 1.5 lc rgb 'black' title "Pressure", \
 "pressure-profile-theoretical_nA00.800000_nB00.800000_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'red' title "Theoretical Pressure"
