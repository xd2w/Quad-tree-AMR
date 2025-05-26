
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


cd $dir

# run the program
time (./FTT_ARM $parFile > out.txt) 2> time.txt 
