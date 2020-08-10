 #!/bin/bash

pushd ..

declare -a dlogz=(
"10.0"
"5.0"
"1.0"
"0.1"
"0.01"
"0.001"
)

declare -a countries=(
"switzerland"
"france"
"germany"
"italy"
)

declare -a models=(
"country.reparam.sird_int.nbin"
"country.reparam.seird_int.nbin"
"country.reparam.seiird2_int.nbin"
"country.reparam.seirud_int.nbin"
"country.reparam.spiird_int.nbin"
)


for model in "${models[@]}"
do
    for c in "${countries[@]}"
    do
        for dlz in "${dlogz[@]}"
        do
            base="./convergence/data/dlogz/$dlz"
            folder="$base/$c/$model"
            mkdir ${folder} -p
 
            msg="first trial, run nested sampling with 1500 samples and dlogz ${dlz}"
            
            outfile="${folder}/knested.out"
            time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
                --silentPlot -ns 1500 -dlz ${dlz} -cm ${model} -c "$c" -ui -ud -df -uint -uip $base -m "${msg}" \
                2>&1 | tee ${outfile}

            python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
    done
done

popd
