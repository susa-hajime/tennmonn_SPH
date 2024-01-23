unset log
reset
set term png
set output "plot.png"
set multiplot

#set size square
set lmargin at screen 0.15
set rmargin at screen 0.97

set xrange[2:8]

set bmargin at screen 0.1
set tmargin at screen 0.54

unset key
set xlabel "x [10^{4}cm]"
set ylabel "rho [10^{-10}g cm^{-3}]"
set yrange[0:0.099]
p  'out.dat' u (($2)*1e-4):($9) w l lw 3, 'out_ini.dat' u 1:2 w l lw 1

set key top left spacing 1.3 Left
set format x ""
set xlabel ""
set ylabel "E_{rad}x10^{14} , e_{gas}x10^{13}  [erg cm^{-3}]"
set bmargin at screen 0.54
set tmargin at screen 0.98
set yrange[0:300]
p  'out.dat' u (($2)*1e-4):($14*1e-13) w l lw 3 title"e_{gas}",'' u (($2)*1e-4):($13*1e-14) w l lw 3 title"E_{rad}", 'out_ini.dat' u 1:($3*1e-14) w l lw 1 notitle, 'out_ini.dat' u 1:($4*1e-13) w l lw 1 notitle
unset multiplot
reset
