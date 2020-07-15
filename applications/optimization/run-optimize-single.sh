 #!/bin/bash

declare -a countries=(
"switzerland"
"france"
"germany"
)

base="./test/"

# for test
model="country.cz_int.nbin"

mkdir ${base} -p

for c in "${countries[@]}"
do
    folder="$base/$c/$model"
    mkdir ${folder} -p

    outfile="${folder}/cmaes.out"
    time PYTHONPATH=../..:../../build:$PYTHONPATH python optimize.py \
        --silentPlot -ns 64 -nt 4 -ng 5000 -cm ${model} -ui -c "$c" -df $base 2>&1 | tee ${outfile}

done
