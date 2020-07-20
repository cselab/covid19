 #!/bin/bash

declare -a countries=(
"china"
#"switzerland"
#"france"
#"germany"
#"italy"
#"russia"
#"us"
#"sweden"
#"brazil"
#"india"
)

base="./data/new/"

#model="country.cz_int.nbin"
model="country.reparam.sird_dint.nbin"

mkdir ${base} -p

for c in "${countries[@]}"
do
    folder="$base/$c/$model"
    mkdir ${folder} -p

    outfile="${folder}/knested.out"
    time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
        --silentPlot -ns 500 -cm ${model} -c "$c" -ui -ud -df $base 2>&1 | tee ${outfile}

    python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
done
