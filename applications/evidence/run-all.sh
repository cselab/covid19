 #!/bin/bash

declare -a arr=(
"country.sir.nbin"
"country.sir.nrm"
"country.sir.tnrm"
"country.sir_int.nbin"
"country.sir_int.nrm"
"country.sir_int.tnrm"
)


for i in "${arr[@]}"
do

   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 2000 -np 2000 -cm "$i"

   folder="data/switzerland/$i/"

   python3 -m korali.plotter --dir "$folder/_korali_samples" --output "$folder/figures/samples.png"

done
