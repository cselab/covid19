 #!/bin/bash

declare -a arr=(
"sir_int"
"seir_int"
"seiir_int"
)

mkdir -p kdata
mkdir -p output_mknested

base="./kdata/"

for model in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot -ns 500 -cm "country.reparam.${model}.tnrm" -c "$c" -df $base --synthetic -dat "./data/${model}_rnd.txt" | tee "./output_mknested/${model}.out"
  
   folder="${base}/country.reparam.${model}.tnrm"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
  
   rm -rf "$folder/_korali_samples"
   rm -rf "$folder/_korali_propagation"
done
