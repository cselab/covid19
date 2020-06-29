 #!/bin/bash

declare -a countries=(
"switzerland"
"france"
"germany"
)

declare -a models=(
"sir_gui"
"seir_gui"
)

base="./gui/nested/"
mkdir ${base} -p

for c in "${countries[@]}"
do
    for m in "${models[@]}"
    do
       model="country.${m}.nbin"
       
       time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
           --silentPlot -ns 1500 -cm ${model} -c ${c} -df $base -nv 45 -pmm False
       
       folder="${base}/${c}/${model}/"

       python3 -m korali.plotter --dir "$folder/_korali_samples" --output "$folder/figures/samples.png"

    done
done
