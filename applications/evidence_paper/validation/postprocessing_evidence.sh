#!/bin/bash


d="./data_0.5_10_1_delay"

outdir="${d}/result/"
mkdir -p $outdir
python3 ./../../../epidemics/utils/postprocessing_nested.py \
    --src "${d}/" --par "D" --out "${outdir}/${m}.csv"

