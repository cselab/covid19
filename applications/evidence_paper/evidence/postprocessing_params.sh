#!/bin/bash

base='/scratch/wadaniel/covid19/intervention/data/g9'
outdir='./figures/'

mkdir -p $outdir

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.sird_int.nbin" "country.reparam.seird_int.nbin" "country.reparam.seiird2_int.nbin" \
    -v "R0" "D" "eps"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.sird_int.nbin" "country.reparam.saphired_int.nbin" "country.reparam.seirud_int.nbin" \
    -v "R0" "D" "eps"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seird_int.nbin" "country.reparam.seiird2_int.nbin" \
    -v "Z"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seirud_int.nbin" "country.reparam.saphired_int.nbin" \
    -v "Zl" "Y" 

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seiird2_int.nbin" "country.reparam.seirud_int.nbin" "country.reparam.saphired_int.nbin" \
    -v "alpha" 

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -c "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -m "country.reparam.seiird2_int.nbin" "country.reparam.saphired_int.nbin" \
    -v "mu" 

exit

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.seiird2_int.nbin" "country.reparam.seirud_int.nbin" "country.reparam.saphire_int.nbin" \
    -v "alpha"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.seirud_int.nbin" "country.reparam.saphire_int.nbin" \
    -v "Zl"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.seirud_int.nbin" "country.reparam.saphire_int.nbin" \
    -v "Y"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.seird_int.nbin" "country.reparam.seiird2_int.nbin" \
    -v "Z"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.sird_int.nbin" "country.reparam.seird_int.nbin" "country.reparam.seiird2_int.nbin" \
    -v "R0" "D"

#python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
#    -m "country.reparam.seirud_int.nbin" "country.reparam.spiird_int.nbin" "country.reparam.seiird2_int.nbin" \
#    -v "alpha"

#python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
#    -m "country.reparam.seirud_int.nbin" "country.reparam.spiird_int.nbin" \
#    -v "Y"

python  ../../../epidemics/utils/plot_comparison.py -df "$base" -sd $outdir \
    -m "country.reparam.sird_int.nbin" "country.reparam.seirud_int.nbin" "country.reparam.saphire_int.nbin" \
    -v "R0" "D"
