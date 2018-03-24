set xlabel "Component A Density"
set ylabel "Component B Density"
set xrange [0.23:0.33]
set yrange [0.23:0.33]
set size square
set key right bottom Left
set key width -10

plot "minimization-unconditionally-2-phase_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' title "Two-phase/Metastable Test Point", \
 "minimization-metastable-densities_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' notitle, \
 "minimization-unconditionally-3-phase_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 13 lc rgb 'dark-grey' title "Three-phase Test Point", \
 "threePhaseDiagram-binodal_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines lc rgb 'black' title "Theoretical Binodal", \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linetype 0 linewidth 2 linecolor rgb 'brown' title "Theoretical Tie Line", \
 "LB-metastable-2-phase-region_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 10 linecolor rgb 'forest-green' title "LB metastable region (2-phase)", \
 "LB-three-phase-region_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 5 linecolor rgb 'black' title "LB 3-phase region", \
 "LB_shield-region-nucleation.dat" with points pointtype 7 linecolor rgb 'yellow' notitle, \
 "singularity-line-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines lc rgb 'black' notitle



