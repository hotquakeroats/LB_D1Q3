set datafile commentschars "#!"
set xlabel "Component A Density"
set ylabel "Component B Density"
set size square
set grid

plot "threePhaseDiagram-densities-twoPhases_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' title "Minimization 2-phase Test Point"
replot "threePhaseDiagram-densities-threePhases_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" using 1:2 with points pointtype 13 linecolor rgb 'dark-grey' title "Minimization 3-phase Test Point"
replot "three-phase-region-coords_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" with points pointtype 7 linecolor rgb 'black' title "3-phase region coordinates"
replot "threePhaseDiagram-tie-lines-theoretical_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" with lines linecolor rgb 'brown' title "Theoretical Tie Line"
replot "threePhaseDiagram-binodal_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" with lines linecolor rgb 'black' title "Theoretical Binodal"
replot "LB-two-phase-region_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" using 1:2 with points pointtype 7 linecolor rgb 'blue' title "LB 2-phase region (log)"
replot "LB-two-phase-region-nid_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" using 1:2 with points pointtype 7 linecolor rgb 'turquoise' title "LB 2-phase region (nid)"
replot "LB-metastable-2-phase-region_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" using 1:2 with points pointtype 7 linecolor rgb 'green' title "LB metastable region (2-phase)"
replot "LB-metastable-3-phase-region_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" using 1:2 with points pointtype 7 linecolor rgb 'red' title "LB metastable region (3-phase)"
replot "LB-three-phase-region_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" using 1:2 with points pointtype 7 linecolor rgb 'gold' title "LB 3-phase region"
replot "singularity-line-coords_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" with lines linecolor rgb 'black' notitle

