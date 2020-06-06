 #!/bin/bash

declare -a arr=(
"sir_raw.txt"
)

base="./data/"

for f in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_nested.py --silentPlot -nt 12 -ns 500 -np 500 -cm "country.sir.tnrm" -c "$c" --sampler 'mTMCMC' -df $base --synthetic -dat $f | tee out-rss.txt
   
   folder="$base/country.sir.tnrm/"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
done
