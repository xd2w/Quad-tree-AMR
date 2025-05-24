
# usage ./run.sh [dir] [parFile]
# if no arguments dir = Run1, parFile = ftt.par

parFile="ftt.par"
dir="Run1"

if [ $# -eq 1 ]; then
    dir=$1
fi

if [ $# -eq 2 ]; then
    dir=$1
    parFile=$2
fi

echo "Running $parFile in $dir"

cd ../
make
cd ./Examples

pwd

mkdir $dir
cd $dir
mkdir DATA

cp ../$parFile ftt.par
cp ../FTT_ARM .

# run the program
# ./FTT_ARM $parFile  
