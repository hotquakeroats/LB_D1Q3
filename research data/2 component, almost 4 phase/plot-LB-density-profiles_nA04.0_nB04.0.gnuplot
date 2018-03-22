set xlabel "Lattice Position" offset 0,0.5,0
set ylabel "Density" offset 2.75,0,0
set key center top

plot  "density-profile-A_0-iterations_nA04.000000_nB04.000000_tcA0.475000_tcB0.475000_aA0.534375_aAB0.301174_aB0.534375.dat" with lines lc rgb 'black' title "Initialized A Density", \
 "density-profile-B_0-iterations_nA04.000000_nB04.000000_tcA0.475000_tcB0.475000_aA0.534375_aAB0.301174_aB0.534375.dat" with lines lc rgb 'red' title "Initialized B Density", \
 "density-profile-A_550000-iterations_nA04.000000_nB04.000000_tcA0.475000_tcB0.475000_aA0.534375_aAB0.301174_aB0.534375.dat" with linespoints pointtype 12 pointsize 1.5 lc rgb 'black' title "Component A Density Profile", \
 "density-profile-B_550000-iterations_nA04.000000_nB04.000000_tcA0.475000_tcB0.475000_aA0.534375_aAB0.301174_aB0.534375.dat" with linespoints pointtype 6 pointsize 1.5 lc rgb 'red' title "Component B Density Profile"
