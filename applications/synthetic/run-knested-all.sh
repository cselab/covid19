 #!/bin/bash

declare -a arr=(
#"sir_int"
#"seir_int"
"seiir_int"
)

#"seir"
#"seiir"
#"seir_int"
#"seiir_int"

#"sir"

mkdir -p kdata
mkdir -p output_mknested

base="./kdata/"

for model in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot -ns 100 -cm "country.reparam.${model}.tnrm" -c "$c" -df $base --synthetic -dat "./data/${model}_rnd.txt" | tee "./output_mknested/${model}.out"
   
done
