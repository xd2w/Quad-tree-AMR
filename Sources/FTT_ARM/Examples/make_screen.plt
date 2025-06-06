set xrange [0:1]
#set yrange [-0.5:0.5]
set yrange [0:1]

set size ratio -1

set style line 1 lw  1

set term eps
set output "out.eps"

plot 'DATA/mesh.000' w l notitle, 'DATA/intf.000' w l lw 5 notitle'