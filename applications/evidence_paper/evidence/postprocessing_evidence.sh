#!/bin/bash

source ../countries.sh

base='/scratch/wadaniel/covid19/intervention/data/g9'
outdir='./result'

mkdir -p $outdir
for c in "${countries[@]}"
do
    python3 ./../../../epidemics/utils/postprocessing_nested.py --src "${base}/${c}/" --par "R0" --out "${outdir}/${c}.csv"

done
