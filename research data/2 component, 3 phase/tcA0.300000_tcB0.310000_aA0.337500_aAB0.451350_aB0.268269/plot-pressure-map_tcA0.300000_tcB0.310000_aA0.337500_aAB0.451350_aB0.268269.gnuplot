set datafile commentschars "#!"
set size square
set grid
stats "pressure-deviations_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" using 3
set cbrange [STATS_max/100:STATS_max]
set logscale cb
set format cb "%.1le%L"

plot "pressure-deviations_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" using 1:2:3 with points pointtype 7 pointsize 0.5 palette title "Pressure Deviations", \
 "threePhaseDiagram-tie-lines-theoretical_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" with lines lc rgb 'black' notitle, \
 "threePhaseDiagram-binodal_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" with lines lc rgb 'black' notitle, \
 "three-phase-region-coords_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" with points pointtype 7 lc rgb 'black' notitle, \
 "singularity-line-coords_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" with lines lc rgb 'black' notitle
