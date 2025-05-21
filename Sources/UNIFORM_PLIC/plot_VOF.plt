set xrange [0:1]
#set yrange [-0.5:0.5]
set yrange [0:1]

set size ratio -1

do for [i=0:9] {
plot  'vof.00'.i u 1:2:3 w image
pause 0.02
}

do for [i=10:99] {
plot  'vof.0'.i u 1:2:3 w image
pause 0.02
}