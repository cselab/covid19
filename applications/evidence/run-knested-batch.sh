 #!/bin/bash

declare -a arr=(
"switzerland"
"france"
)

# OTHER (TOP 10 by Population)

# "russia"
# "germany"
# "turkey"
# "france"
# "united kingdom" @test link
# "italy"
# "spain"
# "ukraine"
# "poland" @test link
# "romania @test link

# "switzerland"
# "sweden"

base="./data/knested/"

declare -a models1=(
"country.reparam.sir_int.nbin"
"country.reparam.sir_int_nogamma.nbin"
"country.reparam.seir_int.nbin"
"country.reparam.seir_int_nogamma.nbin"
"country.reparam.seir_int_nogamma_noZ.nbin"
"country.reparam.seiir_int.nbin" 
"country.reparam.seiir_int_nogamma.nbin" 
"country.reparam.seiir_int_nogamma_noZ.nbin" 
)

declare -a models2=(
"country.reparam.sir_dint.nbin"
"country.reparam.sir_dint_nogamma.nbin"
"country.reparam.seir_dint.nbin"
"country.reparam.seir_dint_nogamma.nbin"
"country.reparam.seir_dint_nogamma_noZ.nbin"
"country.reparam.seiir_dint.nbin" 
"country.reparam.seiir_dint_nogamma.nbin" 
"country.reparam.seiir_dint_nogamma_noZ.nbin" 
)

declare -a models3=(
"country.reparam.sir_int.tnrm"
"country.reparam.sir_int_nogamma.tnrm"
"country.reparam.seir_int.tnrm"
"country.reparam.seir_int_nogamma.tnrm"
"country.reparam.seir_int_nogamma_noZ.tnrm"
"country.reparam.seiir_int.tnrm"
"country.reparam.seiir_int_nogamma.tnrm"
"country.reparam.seiir_int_nogamma_noZ.tnrm"
)

declare -a models4=(
"country.reparam.sir_dint.tnrm"
"country.reparam.sir_dint_nogamma.tnrm"
"country.reparam.seir_dint.tnrm"
"country.reparam.seir_dint_nogamma.tnrm"
"country.reparam.seir_dint_nogamma_noZ.tnrm"
"country.reparam.seiir_dint.tnrm"
"country.reparam.seiir_dint_nogamma.tnrm"
"country.reparam.seiir_dint_nogamma_noZ.tnrm"
)

declare -a models5=(
"country.reparam.sird_int.nbin"
"country.reparam.sird_int_nogamma.nbin"
"country.reparam.seird_int.nbin"
"country.reparam.seird_int_nogamma.nbin"
"country.reparam.seird_int_nogamma_noZ.nbin"
"country.reparam.seiird_int.nbin" 
"country.reparam.seiird_int_nogamma.nbin" 
"country.reparam.seiird_int_nogamma_noZ.nbin" 
)


declare -a models6=(
"country.reparam.sird_dint.nbin"
"country.reparam.sird_dint_nogamma.nbin"
"country.reparam.seird_dint.nbin"
"country.reparam.seird_dint_nogamma.nbin"
"country.reparam.seird_dint_nogamma_noZ.nbin"
"country.reparam.seiird_dint.nbin" 
"country.reparam.seiird_dint_nogamma.nbin" 
"country.reparam.seiird_dint_nogamma_noZ.nbin" 
)


mkdir ${base} -p

for model in "${models6[@]}"
do
    for c in "${arr[@]}"
    do
        folder="$base/$c/$model"
        mkdir ${folder} -p

        outfile="${folder}/knested.out"
        time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
            --silentPlot -ns 500 -cm ${model} -c "$c" -ui -ud -df $base 2>&1 | tee ${outfile}

        python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
done
