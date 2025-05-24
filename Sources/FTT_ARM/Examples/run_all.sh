# for runing all data unattended

./build.sh Lv9_UnifromIntf ftt9.par
./run.sh Lv9_UnifromIntf

echo "Lv9 done"

./build.sh Lv10_UnifromIntf ftt10.par
./run.sh Lv10_UnifromIntf

echo "Lv10 done"

./build.sh Lv11_UnifromIntf ftt11.par
./run.sh Lv11_UnifromIntf

echo "Lv11 done"

# ./build.sh Lv12_5 ftt12.par
# ./run.sh Lv12_5

# echo "Lv12 done"
