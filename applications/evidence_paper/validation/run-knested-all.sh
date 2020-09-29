 #!/bin/bash

declare -a arr=(
#"seirdelay_int"
"seiirdelay_int"
"saphiredelay_int"
"seirudelay_int"
"sirdelay_int"
)

ic=(10)
alpha=(0.3 0.5 0.7)
dispersion=(1)
mkdir -p figures
mkdir -p data

for d in "${dispersion[@]}"
do
for a in "${alpha[@]}"
do
for i0 in "${ic[@]}"
do


base="./data_${a}_${i0}_${d}_delay/"
mkdir -p ${base}

#PYTHONPATH=../../../build:$PYTHONPATH python3 make_synthetic_with_deaths.py \
#    -D 3.7 -Z 3.7 -Zl 0.423 -Y 2.2 -mu 0.4 -I0 ${i0} -alpha ${a} \
#    -dispersion ${d}

PYTHONPATH=../../../build:$PYTHONPATH python3 make_synthetic_with_deaths.py \
    -D 5.2 -Z 5.2 -Zl 0.423 -Y 2.2 -mu 0.4 -I0 ${i0} -alpha ${a} \
    -r0sir 2.25 -r0seir 4.0 -r0seiir 7.0 -r0saphire 1.25 -r0seiru 2.25 -delay 14 \
    -dispersion ${d}

cp ./data/*.txt ${base}
for model in "${arr[@]}"
do

    PYTHONPATH=../../..:../../../build:$PYTHONPATH python sample_knested.py --silentPlot \
       -ns 2500 -dlz 0.1 -bs 12 -nt 12 -cm "country.reparam.${model}.nbin" \
       -ui -ud -uip \
       -c "$c" -df ${base} --synthetic -dat "./${base}/${model}_rnd.txt"
  
   folder="${base}/country.reparam.${model}.nbin"
   python3 -m korali.plotter --dir "$folder/_korali_samples"  --output "$folder/figures/samples.png"
  
done
done
done
done
