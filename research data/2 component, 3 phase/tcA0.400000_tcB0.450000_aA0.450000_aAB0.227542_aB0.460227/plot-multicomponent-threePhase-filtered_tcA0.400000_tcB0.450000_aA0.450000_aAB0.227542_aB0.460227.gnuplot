set datafile commentschars "#!"
set xlabel "Component A Density"
set ylabel "Component B Density"
set size square
plot "minimization-unconditionally-2-phase_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' title "Two-phase/Metastable Test Point", \
 "minimization-metastable-densities_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 13 linecolor rgb 'light-grey' notitle, \
 "minimization-unconditionally-3-phase_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 13 linecolor rgb 'dark-grey' title "Three-phase Test Point", \
 "threePhaseDiagram-binodal_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines linecolor rgb 'black' title "Theoretical Binodal", \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines linetype 0 linewidth 2 linecolor rgb 'brown' title "Theoretical Tie Line", \
 "singularity-line-coords_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines linecolor rgb 'black' notitle, \
 "LB-two-phase-region_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 7 linecolor rgb 'blue' title "LB 2-phase region (log)", \
 "LB-two-phase-region-nid_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 9 linecolor rgb 'light-blue' title "LB 2-phase region (nid)", \
 "LB-metastable-2-phase-region_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 10 linecolor rgb 'web-green' title "LB metastable region (2-phase)", \
 "LB-three-phase-region_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2 with points pointtype 5 linecolor rgb 'black' title "LB 3-phase region", \
