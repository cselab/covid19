 #!/bin/bash

declare -a arr=(
"switzerland"
#"germany"
#"france"
)

base="./data/knested/"
model="country.seir_dint_nogamma_noZ.tnrm"

for c in "${arr[@]}"
do

   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 100000 -np 100000 -cm $model -c "$c" --sampler 'mTMCMC' -df $base | tee -a out-rac.txt

   
   folder="$base/$c/$model"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"

   rm -rf "$folder/_korali_samples"
   rm -rf "$folder/_korali_propagation"

done
