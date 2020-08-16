#!/bin/bash

msg="1 ppm, informed priors, short D"
pushd ..

source countries.sh

name=`whoami`
base="/scratch/${name}/covid19/intervention/data/g9_D52"

declare -a models=(
#"country.reparam.spiird_int.poi"
#"country.reparam.spiird_int.geo"
"country.reparam.spiird_int.nbin"
#"country.reparam.spiird_int.tnrm"
#"country.reparam.spiird_int.tstudent_alt"
)

mkdir ${base} -p

for model in "${models[@]}"
do
    for c in "${countries[@]}"
    do
        folder=$base/$c/$model
        mkdir -p "${folder}"

        outfile=${folder}/knested.out
        time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
            --silentPlot -ns 1500 -dlz 0.1 -cm ${model} -c "$c" -ui -ud -uip -uint -bs 8 -nt 8 -df $base -m "${msg}" \
            2>&1 | tee "${outfile}"

        python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        
        rm -r "$folder/_korali_propagation"

        done
done

popd