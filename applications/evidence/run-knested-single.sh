 #!/bin/bash

msg="test student_t"

declare -a countries=(
#"switzerland"
"france"
#"germany"
#"italy"
#"uk"
#"spain"
#"russia"
#"us"
#"canada"
#"australia"
#"china"
#"japan"
#"south korea"
#"turkey"
)

base="./data/tests/"

model="country.reparam.sird_int.tstudent_alt"
#model="country.cz_int.nbin"
#model="country.reparam.sird_dint.tstudent"
#model="country.reparam.sird_dint.tstudent_alt"
#model="country.reparam.sird_dint.poi"
#model="country.reparam.sird_dint.geo"
#model="country.reparam.sird_dint.tnrm"
#model="country.reparam.sird_dint.nbin"
#model="country.reparam.seird_dint.nbin"
#model="country.reparam.seird_int.nbin"
#model="country.reparam.seiird_dint.nbin"
#model="country.reparam.seiird_int.nbin"
#model="country.reparam.seiird2_int.nbin"
#model="country.reparam.seiird2_dint.nbin"

mkdir ${base} -p

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

