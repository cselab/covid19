 #!/bin/bash

declare -a arr=(
"sir_int_r0_raw"
)

#"sir_int_r0_rndm"


model="country.sir_int_r0_IC.tnrm"

for f in "${arr[@]}"
do
   base="./data/$f/"
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 1000 -np 1000 -cm $model -c "$c" --sampler 'mTMCMC' -df $base --synthetic -dat "$f.txt" | tee out-rsic.txt
   
   folder="$base/$model"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"

   rm -rf "$folder/_korali_samples"
   rm -rf "$folder/_korali_propagation"
done
