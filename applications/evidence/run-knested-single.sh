 #!/bin/bash

msg="test initialization"

declare -a countries=(
#"canada"
#"china"
#"france"
#"germany"
#"italy"
#"japan"
#"russia"
"switzerland"
#"uk"
#"us"
)

base="./data/test_delay/"

#model="country.reparam.seird_int_init.nbin"
#model="country.reparam.seiird2_int_init.nbin"
#model="country.reparam.saphire_int_init.nbin"

#model="country.reparam.sirdelay_int.nbin"
#model="country.reparam.seirdelay_int.nbin"
#model="country.reparam.sird_int.nbin"
#model="country.reparam.seird_int.nbin"
#model="country.reparam.saphired_int.nbin"
#model="country.reparam.saphiredelay_int.nbin"
#model="country.reparam.seirudelay_int.nbin"
model="country.reparam.seirud_int.nbin"
#model="country.reparam.sird_int.nbin"
#model="country.reparam.seiirdelay_int.nbin"
#model="country.reparam.seiird2_int.nbin"

mkdir ${base} -p

for c in "${countries[@]}"
do
    folder=$base/${c}/$model
    mkdir -p "${folder}"

    outfile=${folder}/knested.out
    time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
        --silentPlot -ns 1500 -dlz 0.1 -cm ${model} -c "${c}" -bs 8 -nt 8 -ui -ud -uint -uip -df "$base" -m "${msg}" \
        2>&1 | tee "${outfile}"

    python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
#    rm -rf "$folder/_korali_samples" "$folder/_korali_propagation"
done

