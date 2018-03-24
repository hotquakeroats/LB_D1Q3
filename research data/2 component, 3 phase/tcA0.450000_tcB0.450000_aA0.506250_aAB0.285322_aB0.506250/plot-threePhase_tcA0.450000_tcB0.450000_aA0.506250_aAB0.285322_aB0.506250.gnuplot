set datafile commentschars "#!"
set xlabel "Component A Density"
set ylabel "Component B Density"
set size square
set grid
set key box
plot "threePhaseDiagram-densities-threePhases_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 3 linecolor rgb 'grey' title "3-phase test point"
replot "LB-metastable-3-phase-region_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 3:4 with points pointtype 3 linecolor rgb 'dark-grey' title "Metastable point (3-phase)"
replot "threePhaseDiagram-binodal_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' title "Binodal"
replot "three-phase-region-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'web-green' title "3-phase region"
replot "unconditionally-unstable-region_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'blue' title "Unconditionally unstable"
replot "threePhaseDiagram-densities-threePhases_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 3:4 with points pointtype 1 linecolor rgb 'red' title "Phase 1"
replot "threePhaseDiagram-densities-threePhases_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 5:6 with points pointtype 2 linecolor rgb 'green' title "Phase 2"
replot "threePhaseDiagram-densities-threePhases_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 7:8 with points pointtype 25 linecolor rgb 'black' title "Phase 3"
replot "LB-metastable-3-phase-region_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 7 linecolor rgb 'red' title "LB metastable region (3-phase)"
replot "LB-three-phase-region_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 7 linecolor rgb 'gold' title "LB 3-phase region"
