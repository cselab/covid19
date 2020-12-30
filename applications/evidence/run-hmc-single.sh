 #!/bin/bash

msg="test hmc"

declare -a countries=(
#"canada"
#"china"
"france"
"germany"
#"italy"
#"japan"
#"russia"
"switzerland"
#"uk"
#"us"
)

base="./data/hmc/"

model="country.reparam.sir_int.nbin"
#model="country.reparam.sird_int.nbin"
#model="country.reparam.seird_int_init.nbin"
#model="country.reparam.seiird2_int_init.nbin"
#model="country.reparam.saphire_int_init.nbin"

#model="country.reparam.sirdelay_int.nbin"
#model="country.reparam.seirdelay_int.nbin"
#model="country.reparam.sir_int_nogamma.tnrm"
#model="country.reparam.seird_int.nbin"
#model="country.reparam.saphired_int.nbin"
#model="country.reparam.saphiredelay_int.nbin"
#model="country.reparam.seirudelay_int.nbin"
# model="country.reparam.seirud_int.nbin"
#model="country.reparam.sird_int.nbin"
#model="country.reparam.seiirdelay_int.nbin"
#model="country.reparam.seiird2_int.nbin"

mkdir ${base} -p

for c in "${countries[@]}"
do
    folder=$base/${c}/${model}
    mkdir -p "${folder}"

    outfile=${folder}/knested.out
    time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_hmc.py \
        --silentPlot -cm ${model} -c "${c}" -v "Euclidean" -ns 5000 \
        -ui -uint -df "$base" \
        2>&1 | tee "${outfile}"

    python3 -m korali.plotter --dir "${folder}/_korali_samples"  --output "${folder}/figures/samples.eps" 
    rm -rf "$folder/_korali_samples" "$folder/_korali_propagation"

done

