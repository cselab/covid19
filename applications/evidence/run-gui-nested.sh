 #!/bin/bash

country="france"

declare -a arr=(
"sir_gui"
"seir_gui"
)

base="./gui/nested/"
mkdir ${base} -p

for c in "${arr[@]}"
do
   model="country.${c}.nbin"
   
   time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot -ns 1500 -cm ${model} -c ${country} -df $base
   
   folder="${base}/${country}/${model}/"

   python3 -m korali.plotter --dir "$folder/_korali_samples" --output "$folder/figures/samples.png"

done
