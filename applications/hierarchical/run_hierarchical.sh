 #!/bin/bash


python phase_1.py -m sir_int.nbin -r cantons
python plot_hierarchical.py -m sir_int.nbin -r cantons -ph 1

python phase_2_sir_int.py -m sir_int.nbin -r cantons
python plot_hierarchical.py -m sir_int.nbin -r cantons -ph 2

python phase_3.py -m sir_int.nbin -r cantons
python plot_hierarchical.py -m sir_int.nbin -r cantons -ph 3

python plot_comparison.py -m sir_int.nbin
