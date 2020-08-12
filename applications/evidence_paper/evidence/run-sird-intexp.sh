 #!/bin/bash

msg="run intexp w 1500s, to determine priors"

pushd ..

source countries.sh

name=`whoami`
base="/scratch/${name}/covid19/intervention/data/run1"

declare -a models=(
#"country.reparam.sird_intexp.poi"
"country.reparam.sird_intexp.geo"
"country.reparam.sird_intexp.nbin"
"country.reparam.sird_intexp.tnrm"
#"country.reparam.sird_intexp.tstudent_alt"
)


mkdir ${base} -p

for model in "${models[@]}"
do
    for c in "${countries[@]}"
    do
        folder=$base/${c}/$model
        mkdir -p "${folder}"

        outfile=${folder}/knested.out
        time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
            --silentPlot -ns 1500 -dlz 0.1 -cm ${model} -c "${c}" -bs 8 -nt 8 -ui -ud -df $base -m "${msg}" \
            2>&1 | tee "${outfile}"

        python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
done

popd
