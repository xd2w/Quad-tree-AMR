set xrange [0:1]
#set yrange [-0.5:0.5]
set yrange [0:1]

set size ratio -1

do for [i=0:9] {
plot 'mesh.00'.i w l, 'intf.00'.i w l, 'fgrd.00'.i u 1:2:(-$5*$3*1e-4):(-$5*$4*1e-4):($5) with vector
pause 0.2
}

do for [i=10:99] {
plot 'mesh.0'.i w l, 'intf.0'.i w l, 'fgrd.0'.i u 1:2:(-$5*$3*1e-4):(-$5*$4*1e-4):($5) with vector
pause 0.2
}
