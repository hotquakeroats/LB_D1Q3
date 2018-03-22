reset
unset key
set xlabel "Density"
set ylabel "Temperature"
set zlabel "Free Energy"
set ticslevel 0
set dgrid3d 50,50
set style data lines
set pm3d

splot "free-energy-curves.dat" using 1:2:(-$3)





