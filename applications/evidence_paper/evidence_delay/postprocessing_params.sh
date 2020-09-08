#!/bin/bash

base='/scratch/wadaniel/covid19/data/delay'
outdir='./plots/'

mkdir -p $outdir

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.sirdelay_int.nbin" "country.reparam.seirdelay_int.nbin" "country.reparam.seiirdelay_int.nbin" \
    -v "R0" "D" "eps" "delay"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.sirdelay_int.nbin" "country.reparam.saphiredelay_int.nbin" "country.reparam.seirudelay_int.nbin" \
    -v "R0" "D" "eps" "delay"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seirdelay_int.nbin" "country.reparam.seiirdelay_int.nbin" \
    -v "Z"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seirudelay_int.nbin" "country.reparam.saphiredelay_int.nbin" \
    -v "Zl"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seirudelay_int.nbin" "country.reparam.saphiredelay_int.nbin" \
    -v "Y" 

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seiirdelay_int.nbin" "country.reparam.seirudelay_int.nbin" "country.reparam.saphiredelay_int.nbin" \
    -v "alpha" 

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seiirdelay_int.nbin" "country.reparam.saphiredelay_int.nbin" \
    -v "mu" 
