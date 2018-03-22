set datafile commentschars "#!"
set size square
set grid
stats "chemical-potential-deviations_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 3
set cbrange [STATS_max/1000:STATS_max]
set logscale cb
set format cb "%.1le%L"

plot "chemical-potential-deviations_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" using 1:2:3 with points pointtype 7 pointsize 0.5 palette title "Chemical Potential Deviations - A", \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines lc rgb 'black' notitle, \
 "threePhaseDiagram-binodal_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines lc rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with points pointtype 7 lc rgb 'black' notitle, \
 "singularity-line-coords_tcA0.400000_tcB0.450000_aA0.450000_aAB0.227542_aB0.460227.dat" with lines lc rgb 'black' notitle
