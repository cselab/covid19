 #!/bin/bash

msg="run intlinear w 1500s to determine priors, tact start of intervention, R0 to 30"
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
"turkey"
)

name=`whoami`
base="/scratch/${name}/covid19/intervention/data/run1"

declare -a models=(
#"country.reparam.seird_intsmooth.poi"
"country.reparam.seird_intsmooth.geo"
"country.reparam.seird_intsmooth.nbin"
"country.reparam.seird_intsmooth.tnrm"
#"country.reparam.seird_intsmooth.tstudent_alt"
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
            --silentPlot -ns 1500 -dlz 0.1 -cm ${model} -c "$c" -bs 8 -nt 8 -ui -ud -df $base -m "${msg}" \
            2>&1 | tee "${outfile}"

        python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
done

popd
