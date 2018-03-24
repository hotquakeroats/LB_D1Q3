#set terminal pngcairo size 10in,10in enhanced font "verdana,10"
#set output "LB_shield-region-nucleation.png"
set terminal epslatex size 4.5,4.5 color solid
set output "LB_shield-region-nucleation.eps"

set multiplot layout 3,1

set ylabel "(0.30,0.45)" offset 2.75,0,0
set key outside center above
plot "density-profile-A_1000000-iterations_nA00.300000_nB00.420000_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with linespoints pointtype 12 pointsize 1.5 lc rgb 'black' title "Component A", \
 "density-profile-B_1000000-iterations_nA00.300000_nB00.420000_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with linespoints pointtype 6 pointsize 1.5 lc rgb 'red' title "Component B"

set ylabel "(0.30,1.41)" offset 2.75,0,0
unset key
plot "density-profile-A_1000000-iterations_nA00.300000_nB01.410000_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with linespoints pointtype 12 pointsize 1.5 lc rgb 'black' notitle, \
 "density-profile-B_1000000-iterations_nA00.300000_nB01.410000_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with linespoints pointtype 6 lc rgb 'red' notitle

set xlabel "Lattice Position" offset 0,0.5,0
set ylabel "(0.36,0.52)" offset 2.75,0,0
unset key
plot "density-profile-A_1000000-iterations_nA00.360000_nB00.520000_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with linespoints pointtype 12 pointsize 1.5 lc rgb 'black' notitle, \
 "density-profile-B_1000000-iterations_nA00.360000_nB00.520000_tcA0.450000_tcB0.450000_aA0.506250_aAB0.285322_aB0.506250.dat" with linespoints pointtype 6 pointsize 1.5 lc rgb 'red' notitle

unset multiplot
