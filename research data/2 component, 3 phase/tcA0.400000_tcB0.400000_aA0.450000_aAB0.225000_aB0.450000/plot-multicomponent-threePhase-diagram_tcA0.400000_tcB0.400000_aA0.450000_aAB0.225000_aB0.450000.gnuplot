set datafile commentschars "#!"
set xlabel "Component A Density"
set ylabel "Component B Density"
set size square
set grid

plot "threePhaseDiagram-densities-twoPhases_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' title "Minimization 2-phase Test Point"
replot "threePhaseDiagram-densities-threePhases_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 13 lc rgb 'dark-grey' title "Minimization 3-phase Test Point"
replot "three-phase-region-coords_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with points pointtype 7 lc rgb 'black' title "3-phase region coordinates"
replot "threePhaseDiagram-tie-lines-theoretical_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'brown' title "Theoretical Tie line"
replot "threePhaseDiagram-binodal_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' title "Theoretical Binodal"
replot "LB-two-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 7 lc rgb 'blue' title "LB 2-phase region (log)"
replot "LB-metastable-2-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 7 lc rgb 'green' title "LB metastable region (2-phase)"
replot "LB-metastable-3-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 7 lc rgb 'red' title "LB metastable region (3-phase)"
replot "LB-three-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 7 lc rgb 'gold' title "LB 3-phase region"
replot "singularity-line-coords_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' notitle
