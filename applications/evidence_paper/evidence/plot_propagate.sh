#!/bin/bash

python3 ../../../epidemics/utils/plot_propagation.py \
    --folder "/scratch/wadaniel/covid19/intervention/data/g9" \
    --models "country.reparam.sird_int.nbin" "country.reparam.seird_int.nbin" "country.reparam.saphired_int.nbin" \
    --countries "italy" "switzerland" "uk" "us" "china" \
    -sd "./result_init"


