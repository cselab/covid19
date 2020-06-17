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

base="./data/knested2/"

model="country.reparam.sir_int.tnrm"
# model="country.reparam.seir_int.tnrm"
# model="country.reparam.seiir_int.tnrm"

for c in "${arr[@]}"
do
   time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot -ns 500 -cm ${model} -c "$c" -df $base | tee "knested_${c}_${model}.out"

done
