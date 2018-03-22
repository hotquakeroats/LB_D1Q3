set datafile commentschars "#!"
set xlabel "Component A Density"
set ylabel "Component B Density"
set size square
set grid

plot "LB-metastable-2-phase-region_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 3:4 with points pointtype 13 lc rgb 'light-grey' title "Metastable 2-phase Test Point", \
 "LB-metastable-3-phase-region_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 3:4 with points pointtype 13 lc rgb 'dark-grey' title "Metastable 3-phase Test Point", \
 "three-phase-region-coords_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with points pointtype 7 lc rgb 'black' title "3-phase region coordinates", \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines lc rgb 'brown' title "Theoretical Tie line", \
 "threePhaseDiagram-binodal_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines lc rgb 'black' title "Theoretical Binodal", \
 "LB-metastable-2-phase-region_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 7 lc rgb 'green' title "LB metastable region (2-phase)", \
 "LB-metastable-3-phase-region_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 7 lc rgb 'red' title "LB metastable region (3-phase)", \
 "singularity-line-coords_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines lc rgb 'black' notitle

