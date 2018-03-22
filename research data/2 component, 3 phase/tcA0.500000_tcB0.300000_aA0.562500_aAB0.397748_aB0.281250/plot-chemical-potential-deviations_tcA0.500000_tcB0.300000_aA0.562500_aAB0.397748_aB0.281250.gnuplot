set datafile commentschars "#!"
set key box
plot "chemical-potential-deviations_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 3 with lines title "A deviations"
replot "chemical-potential-deviations_tcA0.500000_tcB0.300000_aA0.562500_aAB0.397748_aB0.281250.dat" using 4 with lines title "B deviations"

