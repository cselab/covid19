 #!/bin/bash

# "switzerland"
# "germany"
declare -a arr=(
"france"
"switzerland"
"germany"
)

# OTHER
# "italy"
# "russia"
# "sweden"
# "ukraine"
# "austria"

base="./data/nested/"
#model="country.sir_int_r0.tnrm"
model="country.reparam.seiir_int.tnrm"


mkdir -p output_nested

for c in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_nested.py --silentPlot -ns 1500 -cm ${model} -c "$c" -df $base

done
