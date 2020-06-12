 #!/bin/bash


declare -a arr=(
"switzerland"
"france"
"germany"
)

# OTHER (TOP 10 by Population)

# "russia"
# "germany"
# "turkey"
# "france"
# "united kingdom" @test link
# "italy"
# "spain"
# "ukraine"
# "poland" @test link
# "romania @test link

# "switzerland"
# "sweden"

base="./data/nested/"

model="country.sir_int_r0.tnrm"
#model="country.seir_int.tnrm"
#model="country.seiir_int.tnrm"


mkdir -p output_nested

for c in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_nested.py --silentPlot -ns 1500 -cm ${model} -c "$c" -df $base
   python plot_nested.py -rf "${base}/${c}/${model}/nested_res.pickle" -of "${base}/${c}/${model}/figures/samples.png"

done
