set datafile commentschars "#!"
set xlabel "Component A Density"
set ylabel "Component B Density"
set size square
set grid
set key box
plot "threePhaseDiagram-densities-threePhases_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 1:2 with points pointtype 3 lc rgb 'grey' title "3-phase test point"
replot "LB-metastable-3-phase-region_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 3:4 with points pointtype 3 lc rgb 'dark-grey' title "Metastable point (3-phase)"
replot "threePhaseDiagram-binodal_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" with lines lc rgb 'black' title "Binodal"
replot "three-phase-region-coords_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" with lines lc rgb 'web-green' title "3-phase region"
replot "unconditionally-unstable-region_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" with lines lc rgb 'blue' title "Unconditionally unstable"
replot "threePhaseDiagram-densities-threePhases_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 3:4 with points pointtype 1 lc rgb 'red' title "Phase 1"
replot "threePhaseDiagram-densities-threePhases_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 5:6 with points pointtype 2 lc rgb 'green' title "Phase 2"
replot "threePhaseDiagram-densities-threePhases_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 7:8 with points pointtype 25 lc rgb 'black' title "Phase 3"
replot "LB-metastable-3-phase-region_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 1:2 with points pointtype 7 lc rgb 'red' title "LB metastable region (3-phase)"
replot "LB-three-phase-region_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 1:2 with points pointtype 7 lc rgb 'gold' title "LB 3-phase region"
