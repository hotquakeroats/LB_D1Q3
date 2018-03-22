set datafile commentschars "#!"
set key box
plot "chemical-potential-deviations_tcA0.350000_tcB0.360000_aA0.393750_aAB0.218725_aB0.337500.dat" using 3 with lines title "A deviations"
replot "chemical-potential-deviations_tcA0.350000_tcB0.360000_aA0.393750_aAB0.218725_aB0.337500.dat" using 4 with lines title "B deviations"

