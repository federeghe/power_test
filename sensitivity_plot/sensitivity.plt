set term png

#
# Various ways to create a 2D heat map from ascii data
#

#set title "Heat Map generated from a file containing Z values only"
#unset key
#set tic scale 0

# Color runs from white to green
#set palette rgbformula -7,2,-7
#set cbrange [0:1]
#set cblabel "Score"
#unset cbtics

set xrange [0:0.5]
set yrange [0:5000]

set xlabel "Shape parameter"
set ylabel "Sample size"

#set multiplot layout 1,2 title "Titolo della figura"
#    margin screen 0.05, 0.95, 0.10, 0.85 spacing screen 0.05

set palette rgbformulae 21, 22, 23

#set view map
#set dgrid3d
#set pm3d interpolate 0,0
#plot "sensitivity0.05.spaces.txt" using 3:2:4 with image pixels

set cbtics ("0" 0, "0.25" 250, "0.5" 500, "0.75" 750, "1" 1000)

unset grid

#set term x11 0

set output "sens-ks-0-01.png"
plot "sensitivity_ks.txt" using 2:1:3 with image pixels

set output "sens-ks-0-05.png"
plot "sensitivity_ks.txt" using 2:1:4 with image pixels

set output "sens-ad-0-01.png"
plot "sensitivity_ad.txt" using 2:1:3 with image pixels

set output "sens-ad-0-05.png"
plot "sensitivity_ad.txt" using 2:1:4 with image pixels

set output "sens-mad-0-01.png"
plot "sensitivity_mad.txt" using 2:1:3 with image pixels

set output "sens-mad-0-05.png"
plot "sensitivity_mad.txt" using 2:1:4 with image pixels

