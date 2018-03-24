set datafile commentschars "#!"
set datafile missing "?"
set size square
set grid

scale = 0.005000
vnorm(x,y) = sqrt(x**2 + y**2)

plot "eigenstuff_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using ($4<=0 ? $1 : 1/0):2:(scale*$6/vnorm($6,1)):(scale/vnorm($6,1)) with vectors linecolor rgb "red" notitle, \
 "eigenstuff_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using ($4>0 ? $1 : 1/0):2:(scale*$6/vnorm($6,1)):(scale/vnorm($6,1)) with vectors linecolor rgb "blue" notitle, \
 "eigenstuff_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using ($5<=0 ? $1 : 1/0):2:(scale*$7/vnorm($7,1)):(scale/vnorm($7,1)) with vectors linecolor rgb "red" notitle, \
 "eigenstuff_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using ($5>0 ? $1 : 1/0):2:(scale*$7/vnorm($7,1)):(scale/vnorm($7,1)) with vectors linecolor rgb "blue" notitle, \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' notitle, \
 "threePhaseDiagram-binodal_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2 with points pointtype 7 linecolor rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 3:4 with points pointtype 7 linecolor rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 5:6 with points pointtype 7 linecolor rgb 'black' notitle, \
 "singularity-line-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' notitle
