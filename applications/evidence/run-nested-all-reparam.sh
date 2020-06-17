 #!/bin/bash


declare -a arr=(
#"switzerland"
#"france"
"germany"
)

# OTHER (TOP 10 by Population)

declare -a arr2=(
"russia"
"turkey"
"france"
"italy"
"spain"
"ukraine"
"poland"
"romania"
"sweden"
"united kingdom"
)

# "germany"

# "switzerland"

base="./data/reparam_intexp/"

# model="country.reparam.sir_int.tnrm"
model="country.reparam.sir_intexp.tnrm"
# model="country.reparam.seir_int.tnrm"
# model="country.reparam.seir_intexp.tnrm"
# model="country.reparam.seiir_int.tnrm"
#model="country.reparam.seiir_intexp.tnrm"


for c in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_nested.py --silentPlot -ns 500 -cm ${model} -c "$c" -df $base
   python plot_nested.py -rf "${base}/${c}/${model}/nested_res.pickle" -of "${base}/${c}/${model}/figures/samples.png"

done
