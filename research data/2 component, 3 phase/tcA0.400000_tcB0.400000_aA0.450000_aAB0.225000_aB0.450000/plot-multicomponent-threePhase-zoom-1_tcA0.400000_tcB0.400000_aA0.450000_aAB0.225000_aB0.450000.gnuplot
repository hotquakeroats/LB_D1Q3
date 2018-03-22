set xlabel "Component A Density"
set ylabel "Component B Density"
set xrange [0.28:1.45]
set yrange [0.28:1.45]
set size square

plot "minimization-unconditionally-2-phase_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' title "Two-phase/Metastable Test Point", \
 "minimization-metastable-densities_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' notitle, \
 "minimization-unconditionally-3-phase_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 13 lc rgb 'dark-grey' title "Three-phase Test Point", \
 "threePhaseDiagram-binodal_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' title "Theoretical Binodal", \
 "unconditionally-unstable-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines linetype 0 linewidth 4 lc rgb 'blue' title "Unconditionally unstable", \
 "three-phase-region-coords_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines linewidth 2 lc rgb 'web-green' title "3-phase region", \
 "singularity-line-coords_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' notitle



