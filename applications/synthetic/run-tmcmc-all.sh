 #!/bin/bash

declare -a arr=(
"sir_int"
#"seir_int"
#"seiir_int"
)

base="./data/tmcmc_rnd_100/"

mkdir -p $base

for model in "${arr[@]}"
do
    PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -ns 5000 -cm "country.reparam.${model}.tstudent_alt" -c "$c" -df $base -nt 8 --synthetic -dat "./data/${model}_rnd.txt" | tee "./${base}/${model}.out"
   
   folder="${base}/country.reparam.${model}.tstudent_alt"
   python3 -m korali.plotter --dir "${folder}/_korali_samples"  --output "${folder}/figures/samples.png"

done
