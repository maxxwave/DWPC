# author: Razvan Ababei

# load the format and color scheme
set term postscript eps enhanced color lw 2.5 font 'Helvetica,18' size 3.5,2.5
#load "~/Desktop/Graphics-tools/ColourPalette.style"

set output  "potential.eps"

set xlabel "X (nm)"
set ylabel "E (10^-18 J)"


p for [i=1:5:1] "f".i."/aici".i u ($3*1e9):($1*1e18) w l ls i  t "".i 

