 #!/bin/bash

declare -a arr=(
#"sird_int"
#"saphired_int"
#"seiird2_int"
"seirud_int"
)

base="./data/"

mkdir -p ${base}

for model in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot \
       -ns 1500 -dlz 0.1 -bs 8 -nt 8 -cm "country.reparam.${model}.nbin" \
       -ui -ud -uip \
       -c "$c" -df $base --synthetic -dat "./data/${model}_rnd.txt"
  
   folder="${base}/country.reparam.${model}.nbin"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
   rm -rf "$folder/_korali_samples" "$folder/_korali_propagation"
  
done
