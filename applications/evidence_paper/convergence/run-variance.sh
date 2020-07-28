 #!/bin/bash

pushd ..

ns=1500
reps=5

msg="first trial, run repeated nested sampling with ${ns} samples, study variance of evidence (${reps} reps)"

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
)


for model in "${models[@]}"
do
    for c in "${countries[@]}"
    do
        for i in {1..5}
        do
            base="./convergence/data/variance/run${i}"
            folder="$base/$c/$model"
            mkdir ${folder} -p
 
            
            outfile="${folder}/knested.out"
            time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
                --silentPlot -ns ${ns} -cm ${model} -dlz 1.0 -bs 8 -nt 8 -c "$c" -ui -ud -df $base -m "${msg}" \
                2>&1 | tee ${outfile}

            python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
    done
done

popd
