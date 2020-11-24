 #!/bin/bash

declare -a arr=(
"sir_dummy"
#"sir_int"
)

base="./data/"

mkdir -p ${base}

for model in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot \
       -ns 1500 -dlz 0.1 -bs 8 -nt 8 -cm "country.reparam.${model}.nbin" \
       -ui \
       -df $base --synthetic -dat "./data/${model}_rnd.txt"
  
    folder="${base}/switzerland/country.reparam.${model}.nbin"
    python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/samples.png"
    rm -rf "$folder/_korali_samples" "$folder/_korali_propagation"
  
done
