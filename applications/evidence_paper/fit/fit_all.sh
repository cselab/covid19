#!/bin/bash

source ../countries.sh

for c in "${countries[@]}"
do
    python3 ../fit_country.py -c "${c}"
done
