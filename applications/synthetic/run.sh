 #!/bin/bash

declare -a arr=(
"sir_int_r0_raw.txt"
"sir_int_r0_rndm.txt"
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

base="./data/"

for f in "${arr[@]}"
do

   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 5000 -np 5000 -cm "country.sir_int_r0.tnrm" -c "$c" --sampler 'mTMCMC' -df $base --synthetic -dat $f
   
   folder="$base/$f/country.sir_int_r0.tnrm/"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"

   rm -rf "$folder/_korali_samples"
   rm -rf "$folder/_korali_propagation"

done
