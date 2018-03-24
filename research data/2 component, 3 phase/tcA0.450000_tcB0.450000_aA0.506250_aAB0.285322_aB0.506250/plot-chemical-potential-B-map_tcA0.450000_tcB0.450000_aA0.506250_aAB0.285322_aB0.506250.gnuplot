set datafile commentschars "#!"
set size square
set grid
stats "chemical-potential-deviations_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 4
set cbrange [STATS_max/1000:STATS_max]
set logscale cb
set format cb "%.1le%L"

plot "chemical-potential-deviations_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 1:2:4 with points pointtype 7 pointsize 0.5 palette title "Chemical Potential Deviations - B"
replot "threePhaseDiagram-tie-lines-theoretical_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' notitle
replot "threePhaseDiagram-binodal_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' notitle
replot "three-phase-region-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with points pointtype 7 linecolor rgb 'black' notitle
replot "singularity-line-coords_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with lines linecolor rgb 'black' notitle