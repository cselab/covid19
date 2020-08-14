 #!/bin/bash

msg="test initialization"

declare -a countries=(
"canada"
"france"
"germany"
"italy"
"switzerland"
)

base="./data/init1/"

declare -a models1=(
"country.reparam.sird_int.nbin"
"country.reparam.seird_int.nbin"
"country.reparam.seirud_int.nbin"
"country.reparam.saphired_int.nbin"
"country.reparam.seiird2_int.nbin"
)


mkdir ${base} -p

for model in "${models1[@]}"
do
    for c in "${countries[@]}"
    do
        folder="$base/$c/$model"
        mkdir ${folder} -p

        outfile="${folder}/knested.out"
        time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
          --silentPlot -ns 1500 -dlz 0.1 -cm ${model} -c "${c}" -bs 8 -nt 8 -ui -ud -uint -uip -df "$base" -m "${msg}" \
        2>&1 | tee "${outfile}"

        python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"

        rm -rf "$folder/_korali_samples"
        rm -rf "$folder/_korali_propagation"

        done
done
