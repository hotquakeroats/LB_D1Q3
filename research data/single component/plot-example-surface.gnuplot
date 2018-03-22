set terminal wxt font "Verdana,14"
unset key
set xlabel "Temperature"
set ylabel "Density"
set zlabel "Free Energy"
set xrange [-4:1]
set yrange [-2:2]
unset xtics
unset ytics
unset ztics
set ticslevel 0
set isosamples 250
set hidden3d
set pm3d
set palette defined (-0.3 "orange", -0.1 "dark-red", -0.02 "blue", 0 "dark-blue")
unset colorbox
set contour base
set view 60,220

f(x,y) = -0.25*(x*x + 2*x*y*y + 2*y**4)

splot f(x,y)





