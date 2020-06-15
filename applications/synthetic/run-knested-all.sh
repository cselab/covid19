 #!/bin/bash

declare -a arr=(
"sir_int_r0"
"seir_int"
"seiir_int"
)

#"seir"
#"seiir"

#"sir"

mkdir -p data
mkdir -p output_knested

base="./data/"

for model in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot -ns 1500 -cm "country.${model}.tnrm" -c "$c" -df $base --synthetic -dat "${model}_rnd.txt" | tee "./output_knested/${model}.out"
   
done
