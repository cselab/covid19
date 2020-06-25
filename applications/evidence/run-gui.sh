 #!/bin/bash

country="france"

declare -a arr=(
"sir_gui"
"seir_gui"
)

base="./gui1/"
mkdir ${base} -p

for c in "${arr[@]}"
do
   model="country.${c}.nbin"
   
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 5000 -np 5000 -c ${country} -cm $model -df ${base} -tc 1.0
   
   folder="${base}/${country}/${model}/"

   python3 -m korali.plotter --dir "$folder/_korali_samples" --output "$folder/figures/samples.png"

done
