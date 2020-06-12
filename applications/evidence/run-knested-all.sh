 #!/bin/bash

#"france"
#"germany"

declare -a arr=(
"switzerland"
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

base="./data/knested/"

model="country.sir_int_r0.tnrm"
#model="country.seir_int.tnrm"
#model="country.seiir_int.tnrm"

for c in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot -ns 3000 -cm ${model} -c "$c" -df $base | tee knested.out
   python plot_nested.py -rf "${base}/${c}/${model}/nested_res.pickle" -of "${base}/${c}/${model}/figures/samples.png"

done
