set datafile commentschars "#!"
set datafile missing "?"
set size square
set grid

scale = 0.005000
vnorm(x,y) = sqrt(x**2 + y**2)

plot "eigenstuff_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using ($4<=0 ? $1 : 1/0):2:(scale*$6/vnorm($6,1)):(scale/vnorm($6,1)) with vectors lc rgb "red" notitle, \
 "eigenstuff_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using ($4>0 ? $1 : 1/0):2:(scale*$6/vnorm($6,1)):(scale/vnorm($6,1)) with vectors lc rgb "blue" notitle, \
 "eigenstuff_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using ($5<=0 ? $1 : 1/0):2:(scale*$7/vnorm($7,1)):(scale/vnorm($7,1)) with vectors lc rgb "red" notitle, \
 "eigenstuff_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using ($5>0 ? $1 : 1/0):2:(scale*$7/vnorm($7,1)):(scale/vnorm($7,1)) with vectors lc rgb "blue" notitle, \
 "threePhaseDiagram-densities-twoPhases_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 1:2 with points pointtype 13 lc rgb 'light-grey' title "Minimization 2-phase Test Point", \
 "threePhaseDiagram-densities-threePhases_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 1:2 with points pointtype 13 lc rgb 'dark-grey' title "Minimization 3-phase Test Point", \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" with lines lc rgb 'black' notitle, \
 "threePhaseDiagram-binodal_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" with lines lc rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 1:2 with points pointtype 7 lc rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 3:4 with points pointtype 7 lc rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 5:6 with points pointtype 7 lc rgb 'black' notitle, \
 "singularity-line-coords_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" with lines lc rgb 'black' notitle
