mkdir -p data
mkdir -p figures
PYTHONPATH=../../build:$PYTHONPATH python3 make_synthetic.py

# for tests
cp data/sir_int_rnd.txt data/sir_dummy_rnd.txt
cp data/sir_int_raw.txt data/sir_dummy_raw.txt
