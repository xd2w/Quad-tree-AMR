
dir="Run1"

if [ $# -eq 1 ]; then
    dir=$1
fi

cd Pastruns
mkdir $dir
cd $dir

cp -r ../DATA .
cp ../ftt.par .