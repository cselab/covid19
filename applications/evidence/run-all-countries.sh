 #!/bin/bash

declare -a arr=(
"switzerland"
"germany"
"france"
"italy"
"russia"
"sweden"
"ukraine"
"austria"
)

# "china"
# "cosovo"

for i in "${arr[@]}"
do

   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 50000 -np 10000 -cm "country.sir_int_r0.tnrm" -c "$i" --sampler 'mTMCMC'
   
   folder="./data/$i/country.sir_int_r0.tnrm/"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"

   rm -rf "$folder/_korali_samples"
   rm -rf "$folder/_korali_propagation"

done
