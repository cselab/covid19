 #!/bin/bash

msg="first trial, run intstep w 500s"
pushd ..

declare -a countries=(
"switzerland"
"france"
)

base="./intervention/data/intstep/"

declare -a models=(
#"country.reparam.sird_ints.poi"
"country.reparam.sird_ints.geo"
"country.reparam.sird_ints.nbin"
"country.reparam.sird_ints.tnrm"
"country.reparam.sird_ints.tstudent"
"country.reparam.sird_ints.tstudent_alt"
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
            --silentPlot -ns 250 -cm ${model} -c "$c" -ui -ud -df $base -m "${msg}" \
            2>&1 | tee ${outfile}

        python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
done

popd
