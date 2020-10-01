#!/bin/bash

base='/scratch/wadaniel/covid19/data/preprocess'
outdir='./params/'

c=("us")
mkdir -p ${outdir}

python3 ./../../../epidemics/utils/postprocessing_params.py --src "${base}/${c}/" --par R0 D Y --out "${outdir}"
