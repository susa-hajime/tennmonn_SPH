unset log
reset
set term png
set output "sub.png"
set multiplot

#set size square
set lmargin at screen 0.15
set rmargin at screen 0.97

set xrange[0.001:2.5]

set bmargin at screen 0.1
set tmargin at screen 0.54

unset key
set xlabel "x [10^{10}cm]"
set ylabel "rho [10^{-10}g cm^{-3}]"
set yrange[0:50]
p  'out.dat' u (($2)*1e-10):($9*1e10) w l lw 3

set key top right spacing 1.3 Left
set format x ""
set xlabel ""
set ylabel "T[K]"
set bmargin at screen 0.54
set tmargin at screen 0.98
set yrange[0:999]
p  'out.dat' u (($2)*1e-10):15 w l lw 3 title"T_{gas}",'' u (($2)*1e-10):16 w l lw 3 title"T_{rad}"
unset multiplot
reset
