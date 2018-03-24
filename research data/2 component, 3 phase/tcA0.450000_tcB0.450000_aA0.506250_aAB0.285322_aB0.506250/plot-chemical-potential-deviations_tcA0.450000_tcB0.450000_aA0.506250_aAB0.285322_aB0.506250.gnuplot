set datafile commentschars "#!"
set key box
plot "chemical-potential-deviations_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 3 with lines title "A deviations"
replot "chemical-potential-deviations_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" using 4 with lines title "B deviations"

