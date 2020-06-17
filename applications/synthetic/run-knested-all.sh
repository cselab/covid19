 #!/bin/bash

declare -a arr=(
"sir_int_r0"
"seir_int"
"seiir_int"
)

#"seir"
#"seiir"
#"seir_int"
#"seiir_int"

#"sir"

mkdir -p data
mkdir -p output_mknested

base="./pdata/"

for model in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot -ns 500 -cm "country.${model}.tnrm" -c "$c" -df $base --synthetic -dat "${model}_rnd.txt" | tee "./output_mknested/${model}.out"
   
done
