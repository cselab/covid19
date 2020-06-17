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

base="./data/knested_intexp/"

# model="country.reparam.sir_intexp.tnrm"
model="country.reparam.seir_intexp.tnrm"
# model="country.reparam.seiir_intexp.tnrm"

for c in "${arr[@]}"
do
   time PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot -ns 500 -cm ${model} -c "$c" -df $base | tee "knested_${c}_${model}.out"

   folder="$base/$c/$model"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
done
