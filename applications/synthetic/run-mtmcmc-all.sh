 #!/bin/bash

declare -a arr=(
"sir_int"
"seir_int"
"seiir_int"
)

base="./mdata/"

mkdir -p $base
mkdir -p output_mtmcmc

for model in "${arr[@]}"
do
    PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -ns 5000 -sa "mTMCMC" -cm "country.reparam.${model}.tnrm" -c "$c" -df $base --synthetic -dat "./data/${model}_rnd.txt" | tee "./output_mtmcmc/${model}.out" -nt 4
   
   folder="${base}/country.reparam.${model}.tnrm"
   python3 -m korali.plotter --dir "${folder}/_korali_samples"  --output "${folder}/figures/samples.png"

done
