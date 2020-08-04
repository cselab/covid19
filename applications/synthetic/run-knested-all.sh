 #!/bin/bash

declare -a arr=(
"sir_int"
"seir_int"
"seiir_int"
)

base="./kdata2/"

for model in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python sample_knested.py --silentPlot -ns 1500 -bs 8 -nt 8 -cm "country.reparam.${model}.tstudent_alt" -c "$c" -df $base --synthetic -dat "./data/${model}_rnd.txt" | tee "${base}/${model}.out" 
  
   folder="${base}/country.reparam.${model}.tstudent_alt"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
  
done
