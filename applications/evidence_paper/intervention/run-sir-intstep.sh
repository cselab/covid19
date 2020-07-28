 #!/bin/bash

msg="run intstep w 1500s to determine priors"
pushd ..

declare -a countries=(
"switzerland"
"france"
"germany"
"italy"
"uk"
"spain"
"russia"
"us"
"canada"
"australia"
"china"
"japan"
"south korea"
)

base="./intervention/data/step/"

declare -a models=(
#"country.reparam.sird_ints.poi"
"country.reparam.sird_ints.geo"
"country.reparam.sird_ints.nbin"
"country.reparam.sird_ints.tnrm"
#"country.reparam.sird_ints.tstudent"
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
            --silentPlot -ns 1500 -dlz 0.1 -cm ${model} -c "$c" -bs 8 -nt 8 -ui -ud -df $base -m "${msg}" \
            2>&1 | tee ${outfile}

        python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
done

popd
