 #!/bin/bash

msg="1 ppm, informed priors, delay, no uint, adjusted alpha prior"
pushd ..

source countries.sh
name=`whoami`

declare -a models=(
"country.reparam.seiirdelay_dint.nbin"
)

for i in {1..1}
do

    base="/scratch/${name}/covid19/data/delay/run_${i}"

    for model in "${models[@]}"
    do
        for c in "${countries[@]}"
        do
            folder=$base/$c/$model
            mkdir -p "${folder}"

            outfile=${folder}/knested.out
            time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py \
                --silentPlot -ns 1500 -dlz 0.1 -cm ${model} -c "$c" -ui -ud -uip -bs 8 -nt 8 -df $base -m "${msg}" \
                2>&1 | tee "${outfile}"

            python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
            
            rm -r "$folder/_korali_propagation"
            
            done
    done

done

popd
