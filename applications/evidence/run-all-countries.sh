 #!/bin/bash

declare -a arr=(
"china"
"germany"
"russia"
"switzerland"
"sweden"
)

for i in "${arr[@]}"
do

   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 25000 -np 5000 -cm "country.sir_int_r0.tnrm" -c "$i" --sampler 'mTMCMC'
   
   folder="./data/$i/country.sir_int_r0.tnrm/"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"

done
