 #!/bin/bash

declare -a arr=( "basic_nbin"
                 "basic_nrm"
                 "basic_tnrm"
                 "altone_nbin"
                 "altone_nrm"
                 "altone_tnrm")


for i in "${arr[@]}"
do

   ./main.py -nt 12 -ns 5000 -cm "sir.$i"

   folder="data/switzerland/sir_$i/"

   python3 -m korali.plotter --dir "$folder/_korali_samples" --output "$folder/figures/samples.png"

done
