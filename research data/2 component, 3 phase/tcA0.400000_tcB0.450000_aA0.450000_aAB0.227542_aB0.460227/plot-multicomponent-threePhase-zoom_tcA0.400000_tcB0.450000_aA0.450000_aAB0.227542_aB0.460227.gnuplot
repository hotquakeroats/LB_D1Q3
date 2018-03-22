set xlabel "Component A Density"
set ylabel "Component B Density"
set xrange [0.7:1.3]
set yrange [0.3:0.7]
set size square
unset key

plot "minimization-unconditionally-2-phase_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' title "Two-phase/Metastable Test Point", \
 "minimization-metastable-densities_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' notitle, \
 "minimization-unconditionally-3-phase_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 13 lc rgb 'dark-grey' title "Three-phase Test Point", \
 "threePhaseDiagram-binodal_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines lc rgb 'black' title "Theoretical Binodal", \
 "LB-metastable-2-phase-region_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 10 linecolor rgb 'web-green' title "LB metastable region (2-phase)", \
 "LB-three-phase-region_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 5 linecolor rgb 'black' title "LB 3-phase region", \
 "singularity-line-coords_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines lc rgb 'black' notitle



