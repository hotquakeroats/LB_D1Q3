set datafile commentschars "#!"
set key box
plot "chemical-potential-deviations_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" using 3 with lines title "A deviations"
replot "chemical-potential-deviations_tcA0.427000_tcB0.427000_aA0.480375_aAB0.270739_aB0.480375.dat" using 4 with lines title "B deviations"

