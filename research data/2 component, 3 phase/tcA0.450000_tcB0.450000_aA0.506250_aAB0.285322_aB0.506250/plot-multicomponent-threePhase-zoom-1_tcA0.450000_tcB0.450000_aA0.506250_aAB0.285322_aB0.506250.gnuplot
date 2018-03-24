set xlabel "Component A Density"
set ylabel "Component B Density"
set xrange [0.2:1.8]
set yrange [0.2:1.8]
set size square


plot "minimization-unconditionally-2-phase_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' title "Two-phase/Metastable Test Point", \
 "minimization-metastable-densities_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' notitle, \
 "minimization-unconditionally-3-phase_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 13 lc rgb 'dark-grey' title "Three-phase Test Point", \
 "threePhaseDiagram-binodal_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines lc rgb 'black' title "Theoretical Binodal", \
 "unconditionally-unstable-region_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linetype 0 linewidth 4 lc rgb 'blue' title "Unconditionally unstable", \
 "three-phase-region-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linewidth 2 lc rgb 'web-green' title "3-phase region", \
 "singularity-line-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines lc rgb 'black' notitle



