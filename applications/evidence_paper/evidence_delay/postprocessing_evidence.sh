#!/bin/bash

source ../countries.sh


for i in {1..1}
do
    base="/scratch/wadaniel/covid19/data/preprocess"
    outdir="./evidence_preprocess"

    mkdir -p $outdir
    for c in "${countries[@]}"
    do
        python3 ./../../../epidemics/utils/postprocessing_nested.py --src "${base}/${c}/" --par "D" --out "${outdir}/${c}_${i}.csv"

    done

done
