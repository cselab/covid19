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

#model="country.seir_int.tnrm"
model="country.seiir_int.tnrm"


mkdir -p output_nested

for c in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_nested.py --silentPlot -ns 1000 -cm ${model} -c "$c" -df $base
   python plot_nested.py -of "${base}/figures/samples.png"

done
