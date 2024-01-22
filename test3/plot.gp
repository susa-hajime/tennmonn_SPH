unset log
reset
set term png
set output "plot.png"
set multiplot

#set size square
set lmargin at screen 0.15
set rmargin at screen 0.97

set xrange[-7:7]

set bmargin at screen 0.1
set tmargin at screen 0.54

unset key
set xlabel "x [10^{14}cm]"
set ylabel "rho [10^{-10}g cm^{-3}]"
set yrange[0.5:3.95]
p  'out_iso.dat' u ($1*1e-14):($2*1e10) w l dt (10,10) title "isothermal", 'out_ad.dat' u ($1*1e-14):($2*1e10) w l dt (10,20) title"adiabatic",'out.dat' u ($2*1e-14):($9*1e10) w l title "data" lw 3


set key top right spacing 1.3 Left
set format x ""
set xlabel ""
set ylabel "T[K]"
set bmargin at screen 0.54
set tmargin at screen 0.98
set yrange[1.05e3:4.5e3]
p  'out_iso.dat' u ($1*1e-14):3 w l dt (10,10) title "isothermal", 'out_ad.dat' u ($1*1e-14):3 w l dt (10,20) title"adiabatic",'out.dat' u ($2*1e-14):15 w l title "data" lw 3

unset multiplot
