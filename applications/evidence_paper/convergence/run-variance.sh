 #!/bin/bash

pushd ..

msg="run repeated nested sampling with 1500 samples, study variance of evidence (5 reps)"

declare -a countries=(
"france"
"germany"
"switzerland"
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
        for i in {1..5}
        do
            base="./convergence/data/variance/run${i}"
            folder="$base/$c/$model"
            mkdir ${folder} -p
 
            
            outfile="${folder}/knested.out"
            time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
                --silentPlot -ns 1500 -cm ${model} -dlz 0.1 -bs 8 -nt 8 -c "$c" -ui -ud -uint -uip -df $base -bs 8 -nt 8 -m "${msg}" \
                2>&1 | tee ${outfile}

            python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
        done
    done
done

popd
