set datafile commentschars "#!"
set size square
set grid
stats "chemical-potential-deviations_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 4
set cbrange [STATS_max/1000:STATS_max]
set logscale cb
set format cb "%.1le%L"

plot "chemical-potential-deviations_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 1:2:4 with points pointtype 7 pointsize 0.5 palette title "Chemical Potential Deviations - B", \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" with lines lc rgb 'black' notitle, \
 "threePhaseDiagram-binodal_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" with lines lc rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" with points pointtype 7 lc rgb 'black' notitle, \
 "singularity-line-coords_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" with lines lc rgb 'black' notitle
