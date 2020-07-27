 #!/bin/bash

msg="first trial, run intlinear w 1500s"
pushd ..

declare -a countries=(
"switzerland"
"france"
"germany"
"italy"
"us"
)

base="./intervention/data/intlinear/"

declare -a models=(
#"country.reparam.sird_int.poi"
"country.reparam.sird_int.geo"
"country.reparam.sird_int.nbin"
"country.reparam.sird_int.tnrm"
"country.reparam.sird_int.tstudent"
"country.reparam.sird_int.tstudent_alt"
)


mkdir ${base} -p

for model in "${models[@]}"
do
    for c in "${countries[@]}"
    do
        folder="$base/$c/$model"
        mkdir ${folder} -p

        outfile="${folder}/knested.out"
        time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
            --silentPlot -ns 1500 -cm ${model} -c "$c" -ui -ud -df $base -m "${msg}" \
            2>&1 | tee ${outfile}

        python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
done

popd
