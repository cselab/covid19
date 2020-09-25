#!/bin/bash

source ../countries.sh


for i in {3..3}
do
    base="/scratch/wadaniel/covid19/data/delay/run_${i}"
    outdir="./result/evidence"

    mkdir -p $outdir
    for c in "${countries[@]}"
    do
        python3 ./../../../epidemics/utils/postprocessing_nested.py --src "${base}/${c}/" --par "D" --out "${outdir}/${c}_${i}.csv"

    done

done
