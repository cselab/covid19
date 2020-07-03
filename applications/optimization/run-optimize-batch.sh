 #!/bin/bash

declare -a arr=(
"switzerland"
#"france"
#"germany"
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

base="./data/optimization/"

# for test
declare -a models0=(
"country.reparam.sir_int.nbin"
)

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
"country.reparam.sir_int.tnrm"
"country.reparam.sir_int_nogamma.tnrm"
"country.reparam.seir_int.tnrm"
"country.reparam.seir_int_nogamma.tnrm"
"country.reparam.seir_int_nogamma_noZ.tnrm"
"country.reparam.seiir_int.tnrm"
"country.reparam.seiir_int_nogamma.tnrm"
"country.reparam.seiir_int_nogamma_noZ.tnrm"
)

declare -a models3=(
"country.reparam.sir_dint.nbin"
"country.reparam.sir_dint_nogamma.nbin"
"country.reparam.seir_dint.nbin"
"country.reparam.seir_dint_nogamma.nbin"
"country.reparam.seir_dint_nogamma_noZ.nbin"
"country.reparam.seiir_dint.nbin" 
"country.reparam.seiir_dint_nogamma.nbin" 
"country.reparam.seiir_dint_nogamma_noZ.nbin" 
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


mkdir ${base} -p

for model in "${models0[@]}"
do
    for c in "${arr[@]}"
    do
        folder="$base/$c/$model"
        mkdir ${folder} -p

        outfile="${folder}/cmaes.out"
        time PYTHONPATH=../..:../../build:$PYTHONPATH python optimize.py \
            --silentPlot -ns 16 -nt 4 -cm ${model} -c "$c" -df $base 2>&1 | tee ${outfile}
        done
done
