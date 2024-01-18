set xrange[0:0.5]
set yrange[0:4.2]
set nokey
set size square
set terminal postscript eps enhanced color size 10cm,10cm
set output "sedov.eps"
p 'out004.dat' u 12:9 with points pt 7 ps 0.2,'analytic/sedov0.04.dat' lw 2 w l
