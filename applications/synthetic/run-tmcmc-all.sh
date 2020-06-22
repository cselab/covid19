 #!/bin/bash

declare -a arr=(
"sir_int"
"seir_int"
"seiir_int"
)

base="./tdata/"

mkdir -p $base
mkdir -p output_tmcmc

for model in "${arr[@]}"
do
    PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -ns 5000 -cm "country.reparam.${model}.tnrm" -c "$c" -df $base -nt 4 --synthetic -dat "./data/${model}_rnd.txt" | tee "./output_tmcmc/${model}.out"
   
   folder="${base}/country.reparam.${model}.tnrm"
   python3 -m korali.plotter --dir "${folder}/_korali_samples"  --output "${folder}/figures/samples.png"

done
