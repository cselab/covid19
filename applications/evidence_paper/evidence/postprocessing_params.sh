#!/bin/bash

base='/scratch/wadaniel/covid19/intervention/data/g9'
outdir='./result/'

mkdir -p $outdir

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.seird_int.nbin" "country.reparam.seiird2_int.nbin" \
    -v "Z"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.sird_int.nbin" "country.reparam.seird_int.nbin" "country.reparam.seiird2_int.nbin" \
    -v "R0" "D" "eps"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.seirud_int.nbin" "country.reparam.spiird_int.nbin" "country.reparam.seiird2_int.nbin" \
    -v "alpha"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.seirud_int.nbin" "country.reparam.spiird_int.nbin" \
    -v "Y"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.sird_int.nbin" "country.reparam.seirud_int.nbin" "country.reparam.spiird_int.nbin" \
    -v "R0" "D" "eps"
