 #!/bin/bash

declare -a arr=(
"sir_int_r0"
"seir"
"seiir"
"seir_int"
"seiir_int"
)

#"sir"

mkdir -p data
mkdir -p output_nested

base="./data/"

for model in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_nested.py --silentPlot -ns 2500 -cm "country.${model}.tnrm" -c "$c" -df $base --synthetic -dat "${model}_raw.txt" | tee "./output_nested/${model}.out"
   
done
