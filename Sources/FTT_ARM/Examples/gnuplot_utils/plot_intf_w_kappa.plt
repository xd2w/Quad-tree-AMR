set size ratio -1
p 'DATA/fgrd.000' u 1:2:($3*($5+1e-50)*1e-5):($4*($5+1e-50)*1e-5) w vector
pause 10