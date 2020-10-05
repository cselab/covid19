#!/bin/bash

python3 ../../../epidemics/utils/plot_propagation.py \
    --folder "/scratch/wadaniel/covid19/data/preprocess" \
    --models "country.reparam.sirdelay_int.nbin" "country.reparam.seirdelay_int.nbin" "country.reparam.seiirdelay_int.nbin" "country.reparam.seirudelay_int.nbin" "country.reparam.saphiredelay_int.nbin" \
    --countries "canada" "china" "france" "germany" "italy" "japan" "russia" "switzerland" "uk" "us" \
    -sd "./propagate_plots"

