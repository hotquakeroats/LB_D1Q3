set xlabel "Component A Density"
set ylabel "Component B Density"
set xrange [1.27:1.43]
set yrange [0.325:0.385]
set size square

plot "threePhaseDiagram-binodal_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' title "Theoretical Binodal", \
 "LB-metastable-2-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 10 linecolor rgb 'web-green' title "LB metastable region (2-phase)", \
 "LB-three-phase-region_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2 with points pointtype 5 linecolor rgb 'black' title "LB 3-phase region", \
 "singularity-line-coords_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' notitle



