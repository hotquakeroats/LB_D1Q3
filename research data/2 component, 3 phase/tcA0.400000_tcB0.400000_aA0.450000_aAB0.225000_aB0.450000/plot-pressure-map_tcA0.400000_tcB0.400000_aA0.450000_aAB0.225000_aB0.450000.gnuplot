set datafile commentschars "#!"
set size square
set grid
stats "pressure-deviations_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 3
set cbrange [STATS_max/100:STATS_max]
set logscale cb
set format cb "%.1le%L"

plot "pressure-deviations_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" using 1:2:3 with points pointtype 7 pointsize 0.5 palette title "Pressure Deviations"
replot "threePhaseDiagram-tie-lines-theoretical_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' notitle
replot "threePhaseDiagram-binodal_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' notitle
replot "three-phase-region-coords_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with points pointtype 7 lc rgb 'black' notitle
replot "singularity-line-coords_tcA0.400000_tcB0.400000_aA0.450000_aAB0.225000_aB0.450000.dat" with lines lc rgb 'black' notitle
