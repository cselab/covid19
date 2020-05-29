 #!/bin/bash

declare -a arr=(
"sir_raw.txt"
)


base="./data/"

for f in "${arr[@]}"
do
   PYTHONPATH=../..:../../build:$PYTHONPATH python optimize.py --silentPlot -nt 12 -ns 16 -mg 2000 -cm "country.sir.tnrm" -c "$c" -df $base --synthetic -dat $f
   
done
