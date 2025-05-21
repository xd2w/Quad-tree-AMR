set xrange [0:1]
#set yrange [-0.5:0.5]
set yrange [0:1]

set size ratio -1

do for [i=0:9] {
plot 'intf.00'.i w l
pause 0.2
}

do for [i=10:99] {
plot 'intf.0'.i w l
pause 0.2
}
