 #!/bin/bash

declare -a arr=(
"sir_int_r0_rndm"
"sir_int_r0_raw"
)



model="country.sir_int_r0.tnrm"

for f in "${arr[@]}"
do
   base="./data/$f/"
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 25000 -np 25000 -cm $model -c "$c" --sampler 'mTMCMC' -df $base --synthetic -dat "$f.txt" | tee out-rsi.txt
   
   folder="$base/$model"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"

   rm -rf "$folder/_korali_samples"
   rm -rf "$folder/_korali_propagation"
done
