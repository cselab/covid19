mkdir -p data
mkdir -p figures
PYTHONPATH=../../build:$PYTHONPATH python3 make_synthetic_with_deaths.py
