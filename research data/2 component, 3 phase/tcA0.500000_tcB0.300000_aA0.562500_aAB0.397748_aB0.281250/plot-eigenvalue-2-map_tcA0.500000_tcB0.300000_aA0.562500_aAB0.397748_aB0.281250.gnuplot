set datafile commentschars "#!"
set datafile missing "?"
set size square
set grid
stats "eigenstuff_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 5
set cbrange [-1*log(abs(STATS_min)+1):log(abs(STATS_max)+1)]
set palette defined (-1*log(abs(STATS_min)+1) "orange", -1e-6 "red", 0 "white", 1e-6 "blue", log(abs(STATS_max)+1) "black")

plot "eigenstuff_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 1:2:($5<0 ? -1*log(abs($5)+1) : log(abs($5)+1)) with points pointtype 7 pointsize 1.5 palette title "Eigenvalue 2", \
 "threePhaseDiagram-densities-twoPhases_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' title "Minimization 2-phase Test Point", \
 "threePhaseDiagram-densities-threePhases_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 1:2 with points pointtype 13 lc rgb 'dark-grey' title "Minimization 3-phase Test Point", \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" with lines lc rgb 'black' notitle, \
 "threePhaseDiagram-binodal_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" with lines lc rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 1:2 with points pointtype 7 lc rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 3:4 with points pointtype 7 lc rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 5:6 with points pointtype 7 lc rgb 'black' notitle, \
 "singularity-line-coords_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" with lines lc rgb 'black' notitle
