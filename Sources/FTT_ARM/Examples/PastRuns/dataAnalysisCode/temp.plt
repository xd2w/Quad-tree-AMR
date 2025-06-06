
set term eps
set output "temp_out_intf.400.eps"
        
set size ratio -1

set yrange [0.7:0.9]
set xrange [0:0.2]
p '../Lv12_3HF/DATA/thintf.200' w l notitle, '../Lv10_8HF/DATA/intf.200' w l notitle, 
