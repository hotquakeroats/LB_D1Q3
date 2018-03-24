set datafile commentschars "#!"
set key box
plot "chemical-potential-deviations_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 3 with lines title "A deviations"
replot "chemical-potential-deviations_tcA0.320000_tcB0.320000_aA0.360000_aAB0.072000_aB0.360000.dat" using 4 with lines title "B deviations"

