#!/bin/bash

base='/scratch/wadaniel/covid19/data/preprocess'
outdir='./comparison_preprocess/'

mkdir -p ${outdir}

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "switzerland" "uk" "us" \
    -m "country.reparam.sirdelay_int.nbin" "country.reparam.seirdelay_int.nbin" "country.reparam.seiirdelay_int.nbin" "country.reparam.seirudelay_int.nbin" "country.reparam.saphiredelay_int.nbin" \
    -v "r" "delay" "tact" "dtact" "kbeta" "R0" "D" "eps" 

exit

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.sirdelay_int.nbin" "country.reparam.seirdelay_int.nbin" "country.reparam.seiirdelay_int.nbin" "country.reparam.seirudelay_int.nbin" "country.reparam.saphiredelay_int.nbin" \ 
    -v "r" "delay" "tact" "dtact" "kbeta" "R0" "D" "eps" 

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seirdelay_int.nbin" "country.reparam.seiirdelay_int.nbin" \
    -v "Z"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seirudelay_int.nbin" "country.reparam.saphiredelay_int.nbin" \
    -v "Zl" "Y"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seiirdelay_int.nbin" "country.reparam.seirudelay_int.nbin" "country.reparam.saphiredelay_int.nbin" \
    -v "alpha" "mu"
