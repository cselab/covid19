 #!/bin/bash

declare -a arr=(
"switzerland"
"sweden"
"germany"
"russia"
"china"
)


for i in "${arr[@]}"
do

   PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py --silentPlot -nt 12 -ns 25000 -np 5000 -cm "country.sir_int_r0.tnrm" -c "$i" --sampler 'mTMCMC'

done
