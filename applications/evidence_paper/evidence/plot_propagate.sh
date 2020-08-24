#!/bin/bash

python3 ../../../epidemics/utils/plot_propagation.py \
    --folder "/scratch/wadaniel/covid19/intervention/data/g9" \
    --models "country.reparam.sird_int.nbin" "country.reparam.seird_int.nbin" \
    --countries "uk" "italy" "switzerland" \
    -sd "./result"

exit

python3 ../../../epidemics/utils/plot_propagation.py \
    --folder "/scratch/wadaniel/covid19/intervention/data/g9" \
    --models "country.reparam.sird_int.nbin" "country.reparam.seiird2_int.nbin" \
    --countries "uk" "italy" "switzerland" \
    -sd "./result"

python3 ../../../epidemics/utils/plot_propagation.py \
    --folder "/scratch/wadaniel/covid19/intervention/data/g9" \
    --models "country.reparam.sird_int.nbin" "country.reparam.saphired_int.nbin" \
    --countries "uk" "italy" "switzerland" \
    -sd "./result"

python3 ../../../epidemics/utils/plot_propagation.py \
    --folder "/scratch/wadaniel/covid19/intervention/data/g9" \
    --models "country.reparam.sird_int.nbin" "country.reparam.seirud_int.nbin" \
    --countries "uk" "italy" "switzerland" \
    -sd "./result"
