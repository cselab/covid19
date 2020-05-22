#PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py -nt 12 -ns 5000 -cm "country.sir_int_r0.tnrm" -c "china" --sampler 'mTMCMC'
#PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py -nt 12 -ns 5000 -cm "country.sir_int_r0.tnrm" -c "germany" --sampler 'mTMCMC'
#PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py -nt 12 -ns 5000 -cm "country.sir_int_r0.tnrm" -c "russia" --sampler 'mTMCMC'
#PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py -nt 12 -ns 5000 -cm "country.sir_int_r0.tnrm" -c "sweden" --sampler 'mTMCMC'
#PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py -nt 12 -ns 5000 -cm "country.sir_int.tnrm" -c "switzerland" --sampler 'mTMCMC'
#PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py -nt 12 -ns 5000 -cm "country.sir_gui.tnrm" -c "switzerland" --sampler 'mTMCMC'
PYTHONPATH=../..:../../build:$PYTHONPATH python sample.py -nt 12 -ns 100 -cm "country.sir.nbin" -c "switzerland"
