 #!/bin/bash

declare -a arr=(
"sir.nbin"
"sir.nrm"
"sir.tnrm"
"sir_int.nbin"
"sir_int.nrm"
"sir_int.tnrm"
)


for i in "${arr[@]}"
do

   ./sample.py --silentPlot -nt 12 -ns 2000 -cm "$i"

   folder="data/switzerland/$i/"

   python3 -m korali.plotter --dir "$folder/_korali_samples" --output "$folder/figures/samples.png"

done
