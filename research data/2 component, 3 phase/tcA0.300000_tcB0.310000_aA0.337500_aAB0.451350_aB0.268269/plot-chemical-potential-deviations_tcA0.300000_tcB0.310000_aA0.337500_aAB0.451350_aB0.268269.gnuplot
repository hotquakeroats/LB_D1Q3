set datafile commentschars "#!"
set key box
plot "chemical-potential-deviations_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" using 3 with lines title "A deviations"
replot "chemical-potential-deviations_tcA0.300000_tcB0.310000_aA0.337500_aAB0.451350_aB0.268269.dat" using 4 with lines title "B deviations"

