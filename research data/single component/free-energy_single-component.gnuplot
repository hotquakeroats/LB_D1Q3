set terminal wxt font "Verdana,14"
set xlabel "Density"
set ylabel "Free Energy Density"
set xtics ("binodal" 0, "spinodal" 0.2, "spinodal" 0.81, "binodal" 1)
unset ytics
unset key

set label "Stable" at -0.18,-1.15 center
set object 1 rect from graph 0.2168674699, graph 0 to graph 0.3373493976, graph 1 back
set object 1 rect fc rgb "light-gray" fillstyle solid 0.15 noborder
set label "Metatable" at 0.1,-1.15 center
set object 2 rect from graph 0.3373493976, graph 0 to graph 0.7048192771, graph 1 back
set object 2 rect fc rgb "light-gray" fillstyle solid 0.25 noborder
set label "Unconditionally Unstable" at 0.5,-1.15 center
set object 3 rect from graph 0.7048192771, graph 0 to graph 0.8192771084, graph 1 back
set object 3 rect fc rgb "light-gray" fillstyle solid 0.1 noborder
set label "Metatable" at 0.9,-1.15 center
set label "Stable" at 1.15,-1.15 center

f(x) = a*x + b
fit f(x) "free-energy-curve_binodal-points.dat" using 3:4 via a,b

plot f(x) lc rgb 'black', \
 "free-energy-curve_stable-regions.dat" using 3:4 with lines lc rgb 'blue', \
 "free-energy-curve_metastable-regions.dat" using 3:4 with lines lc rgb 'green', \
 "free-energy-curve_unstable-region.dat" using 3:4 with lines lc rgb 'red', \
 "free-energy-curve_spinodal-points.dat" using 3:4 with points pointtype 7 lc rgb 'black', \
 "free-energy-curve_binodal-points.dat" using 3:4 with points pointtype 13 lc rgb 'black'

