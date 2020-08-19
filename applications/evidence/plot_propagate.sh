#!/bin/bash

python3 ../../epidemics/utils/plot_propagation.py --folder "./data/test/" \
    --models "country.reparam.sird_int.nbin" "country.reparam.seird_int.nbin" "country.reparam.saphired_int.nbin" \
    --countries "canada" "france" -sd "."


