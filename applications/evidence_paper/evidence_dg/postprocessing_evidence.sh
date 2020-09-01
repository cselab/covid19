#!/bin/bash

source ../countries.sh

base='/scratch/wadaniel/covid19/intervention/data/g9_new'
outdir='./result_g9_new'

mkdir -p $outdir
for c in "${countries[@]}"
do
    python3 ./../../../epidemics/utils/postprocessing_nested.py --src "${base}/${c}/" --par "D" --out "${outdir}/${c}.csv"

done
