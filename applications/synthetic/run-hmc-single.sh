 #!/bin/bash

msg="test hmc with synthetic data"

declare -a arr=(
"sir_int"
#"sir_dummy"
)

base="./data/"

mkdir -p ${base}

for model in "${arr[@]}"
do
    PYTHONPATH=../..:../../build:$PYTHONPATH \
        python sample_hmc.py \
        --silentPlot -cm ${model} -cm "country.reparam.${model}.nbin" -v "Euclidean" -ns 1500 \
        -ui \
        -df "$base" --synthetic -dat "./data/${model}_rnd.txt"
   
    folder="${base}/switzerland/country.reparam.${model}.nbin"
    python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/samples.png"
    rm -rf "$folder/_korali_samples" "$folder/_korali_propagation"
done

