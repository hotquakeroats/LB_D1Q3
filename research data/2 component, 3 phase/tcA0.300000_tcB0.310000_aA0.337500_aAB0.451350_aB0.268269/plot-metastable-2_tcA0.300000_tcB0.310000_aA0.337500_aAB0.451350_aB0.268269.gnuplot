set datafile commentschars "#!"
set xlabel "Component A Density"
set ylabel "Component B Density"
set size square
set grid

plot "LB-metastable-2-phase-region_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" using 3:4 with points pointtype 13 lc rgb 'light-grey' title "Metastable 2-phase Test Point", \
 "three-phase-region-coords_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" with points pointtype 7 lc rgb 'black' title "3-phase region coordinates", \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" with lines lc rgb 'brown' title "Theoretical Tie line", \
 "threePhaseDiagram-binodal_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" with lines lc rgb 'black' title "Theoretical Binodal", \
 "LB-metastable-2-phase-region_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" using 1:2 with points pointtype 7 lc rgb 'green' title "LB metastable region (2-phase)", \
 "singularity-line-coords_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" with lines lc rgb 'black' notitle

