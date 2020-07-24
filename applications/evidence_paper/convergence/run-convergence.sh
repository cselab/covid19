 #!/bin/bash

pushd ..

declare -a nsamples=(
"100"
"500"
"1000"
"1500"
"2000"
"3000"
"5000"
)

declare -a countries=(
"switzerland"
"france"
"germany"
"italy"
)

declare -a models=(
"country.reparam.sird_intsmooth.nbin"
"country.reparam.seird_intsmooth.nbin"
"country.reparam.seiird2_intsmooth.nbin"
)


mkdir ${base} -p

   
for model in "${models[@]}"
do
    for c in "${countries[@]}"
    do
        for ns in "${nsamples[@]}"
        do
            base="./convergence/data/$ns"
            folder="$base/$c/$model"
            mkdir ${folder} -p
 
            msg="first trial, run nested sampling with ${ns} samples"
            
            outfile="${folder}/knested.out"
            time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
                --silentPlot -ns ${ns} -cm ${model} -c "$c" -ui -ud -df $base -m "${msg}" \
                2>&1 | tee ${outfile}

            python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
    done
done

popd
