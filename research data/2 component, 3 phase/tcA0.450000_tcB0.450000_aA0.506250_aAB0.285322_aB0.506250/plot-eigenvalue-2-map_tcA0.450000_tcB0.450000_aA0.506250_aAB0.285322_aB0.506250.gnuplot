set datafile commentschars "#!"
set datafile missing "?"
set size square
set grid
stats "eigenstuff_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 5
set cbrange [-1*log(abs(STATS_min)+1):log(abs(STATS_max)+1)]
set palette defined (-1*log(abs(STATS_min)+1) "orange", -1e-6 "red", 0 "white", 1e-6 "blue", log(abs(STATS_max)+1) "black")

plot "eigenstuff_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2:($5<0 ? -1*log(abs($5)+1) : log(abs($5)+1)) with points pointtype 7 pointsize 1.5 palette title "Eigenvalue 2"
replot "threePhaseDiagram-densities-twoPhases_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' title "Minimization 2-phase Test Point"
replot "threePhaseDiagram-densities-threePhases_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 13 linecolor rgb 'dark-grey' title "Minimization 3-phase Test Point"
replot "threePhaseDiagram-tie-lines-theoretical_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' notitle
replot "threePhaseDiagram-binodal_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' notitle
replot "three-phase-region-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 7 linecolor rgb 'black' notitle
replot "three-phase-region-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 3:4 with points pointtype 7 linecolor rgb 'black' notitle
replot "three-phase-region-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 5:6 with points pointtype 7 linecolor rgb 'black' notitle
replot "singularity-line-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' notitle
