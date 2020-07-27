 #!/bin/bash

pushd ..

ns=1500

declare -a batch=(
"1"
"2"
"4"
"8"
"16"
"32"
)

declare -a countries=(
"france"
"germany"
"italy"
)

declare -a models=(
"country.reparam.sird_int.nbin"
"country.reparam.seird_int.nbin"
"country.reparam.seiird2_int.nbin"
)


for model in "${models[@]}"
do
    for c in "${countries[@]}"
    do
        for b in "${batch[@]}"
        do
            base="./convergence/data/strong_batch/$b"
            folder="$base/$c/$model"
            mkdir ${folder} -p
 
            msg="first trial, run nested sampling with ${ns} batch size 1500 samples, 0.1 dlogz"
            
            outfile="${folder}/knested.out"
            time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
                --silentPlot -ns ${ns} -bs ${b} -cm ${model} -c "$c" -nt 12 -ui -ud -df $base -m "${msg}" \
                2>&1 | tee ${outfile}

            python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
    done
done

popd
