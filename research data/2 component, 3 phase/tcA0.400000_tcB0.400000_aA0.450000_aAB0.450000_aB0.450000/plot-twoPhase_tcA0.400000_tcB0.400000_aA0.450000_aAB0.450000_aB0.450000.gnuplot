set datafile commentschars "#!"
set xlabel "Component A Density"
set ylabel "Component B Density"
set size square
set grid
set key box
plot "threePhaseDiagram-densities-twoPhases_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" using 1:2 with points pointtype 3 lc rgb 'grey' title "2-phase test point"
replot "LB-metastable-2-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" using 3:4 with points pointtype 3 lc rgb 'dark-grey' title "Metastable point (2-phase)"
replot "LB-metastable-3-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" using 3:4 with points pointtype 3 lc rgb 'white' notitle
replot "threePhaseDiagram-tie-lines-theoretical_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" with lines lc rgb 'red' title "Tie line"
replot "threePhaseDiagram-binodal_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" with lines lc rgb 'black' title "Binodal"
replot "three-phase-region-coords_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" with lines lc rgb 'web-green' title "3-phase region"
replot "threePhaseDiagram-densities-twoPhases_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" using 3:4 with points pointtype 1 lc rgb 'red' title "Phase 1"
replot "threePhaseDiagram-densities-twoPhases_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" using 5:6 with points pointtype 2 lc rgb 'black' title "Phase 2"
replot "LB-metastable-2-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" using 1:2 with points pointtype 7 lc rgb 'green' title "LB metastable region (2-phase)"
replot "LB-two-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" using 1:2 with points pointtype 7 lc rgb 'blue' title "LB 2-phase region"
