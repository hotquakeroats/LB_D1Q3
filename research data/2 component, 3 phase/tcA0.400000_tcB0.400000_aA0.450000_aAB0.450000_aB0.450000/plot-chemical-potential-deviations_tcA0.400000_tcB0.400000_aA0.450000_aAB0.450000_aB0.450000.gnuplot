set datafile commentschars "#!"
set key box
plot "chemical-potential-deviations_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" using 3 with lines title "A deviations"
replot "chemical-potential-deviations_tcA0.400000_tcB0.400000_aA0.450000_aAB0.450000_aB0.450000.dat" using 4 with lines title "B deviations"

