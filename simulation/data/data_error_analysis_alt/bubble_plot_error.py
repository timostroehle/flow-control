### f
# set terminal epslatex lw 2 size 8cm,7cm "Helvetica,8pt" #color colortext
# set output './figures/error.tex'
 
set xlabel '\(n\)' #offset 0,1,0
set ylabel '\(\varepsilon\)' #offset 3,0,0

set border lw 1
set key at graph 1.0,0.60 spacing 1.5

set yrange [1e-6:3e-1]
set format y "10^{%.0T}"
set ytics 1e-6,1e-1 logscale

set xrange [1.8:14]
set format x "%.0f"
set xtics (0,2,4,6,8,10,12) 


set logscale xy
set grid ytics, xtics, mxtics, mytics

set arrow 1 from 2e0,1e-5 to 4e0,1e-5 nohead
set arrow 2 from 4e0,1e-5 to 2e0,4e-5 nohead
set arrow 3 from 2e0,4e-5 to 2e0,1e-5 nohead 
set label "2" at 2.1e0,1.7e-5

set arrow 4 from 2e0,1e-2 to 4e0,1e-2 nohead
set arrow 5 from 4e0,1e-2 to 2e0,2e-2 nohead
set arrow 6 from 2e0,2e-2 to 2e0,1e-2 nohead 
set label "1" at 2.1e0,1.35e-2

plot './data/data_error_analysis/err_mid.dat' using 1:2 w lp lt 8 pt 6 ps 2 lw 1 dt 1 title 'Mid-point',\
     './data/data_error_analysis/err_imp.dat' using 1:2 w lp lt 8 pt 8 ps 2.5 lw 1 dt 1 title 'Euler impl.'
     
pause -1



