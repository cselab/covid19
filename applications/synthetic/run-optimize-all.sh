 #!/bin/bash

declare -a arr=(
"sir"
"sir_int_r0"
"seir"
"seir_int"
"seiir"
"seiir_int"
)

mkdir -p data
mkdir -p output

base="./data/"

for model in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python optimize.py --silentPlot -nt 12 -ns 32 -mg 5000 -cm "country.${model}.tnrm" -c "$c" -df $base --synthetic -dat "./data/${model}_raw.txt" | tee "./output/${model}.out"
   
done
