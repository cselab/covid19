#!/bin/bash

base='/scratch/wadaniel/covid19/data/preprocess'

source ../countries.sh

for c in "${countries[@]}"
do  
    outdir="./params/${c}"
    mkdir -p ${outdir}
    python3 ./../../../epidemics/utils/postprocessing_params.py --src "${base}/${c}/" --par R0 D Z Zl Y alpha mu --out "${outdir}"
done
