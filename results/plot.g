# author: Razvan Ababei

# load the format and color scheme
set term postscript eps enhanced color lw 1.5 font 'Helvetica,18' size 3.5,2.5
#load "~/Desktop/Graphics-tools/ColourPalette.style"

set output  "Bifurcation.eps"
set multiplot layout 2,2 columnsfirst scale 1.1,0.9
unset xlabel 
unset ylabel
unset xtics
unset ytics

p "f1/Bifurcation.data" u 1:($2*1e9) w p ps 0.2 pt 1 t "1"

p "f2/Bifurcation.data" u 1:($2*1e9) w p ps 0.2 pt 2 t "2"

p "f3/Bifurcation.data" u 1:($2*1e9) w p ps 0.2 pt 3 t "3"

p "f4/Bifurcation.data" u 1:($2*1e9) w p ps 0.2 pt 4 t "4"

unset multiplot
