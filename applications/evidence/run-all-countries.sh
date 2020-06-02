 #!/bin/bash

declare -a arr=(
"switzerland"
"germany"
"france"
)

# OTHER
# "italy"
# "russia"
# "sweden"
# "ukraine"
# "austria"

# NOT WORKING
# "china"
# "cosovo"

base="./data/nogamma/"
model="country.sir_int_r0_nogamma.tnrm"

for c in "${arr[@]}"
do

   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 25000 -np 25000 -cm $model -c "$c" --sampler 'mTMCMC' -df $base

   
   folder="$base/$c/$model"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"

   rm -rf "$folder/_korali_samples"
   rm -rf "$folder/_korali_propagation"

done
