#!/bin/bash

source ../countries.sh


for i in {1..1}
do
    base="/scratch/wadaniel/covid19/data/stat/run_${i}"
    outdir="./result/stat/evidence"

    mkdir -p $outdir
    for c in "${countries[@]}"
    do
        python3 ./../../../epidemics/utils/postprocessing_nested.py --src "${base}/${c}/" --par "D" --out "${outdir}/${c}_${i}.csv"

    done

done
