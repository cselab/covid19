#!/bin/bash

source ../countries.sh

base='/scratch/wadaniel/covid19/data/delay'
outdir='./result/delay'

mkdir -p $outdir
for c in "${countries[@]}"
do
    python3 ./../../../epidemics/utils/postprocessing_nested.py --src "${base}/${c}/" --par "D" --out "${outdir}/${c}.csv"

done
