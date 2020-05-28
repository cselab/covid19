 #!/bin/bash

declare -a arr=(
"sir_int_r0_raw.txt"
)

#"sir_int_r0_rndm.txt"


base="./data/"

for f in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python optimize.py --silentPlot -nt 12 -ns 32 -mg 1000 -cm "country.sir_int_r0.tnrm" -c "$c" -df $base --synthetic -dat $f
   
#   folder="$base/country.sir_int_r0.tnrm/"
#   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
done
