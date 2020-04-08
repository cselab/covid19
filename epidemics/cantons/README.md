# Multi-regional SEII model.

Reference:

Substantial undocumented infection facilitates the rapid dissemination of novel coronavirus (SARS-CoV2), Li et al (2020), https://science.sciencemag.org/content/early/2020/03/24/science.abb3221/

The model keeps track of 5 values for each canton:
- S -susceptible
- E - exposed
- Ir - documented (registered) infected
- Iu - undocumented infected
- N - canton population

The paper does not treat N as a constant.
However, it is reasonable to assume that migrations are negligible in the time span of an infection.
Thus, we symmetrize the matrix `Mij` by modifying it as `M'ij = Mij + Mji`, which accounts for the fact that people living in `i` and working in `j` also have to go back.

Internally, the state is stored as a vector of 5 * 26 values, `[S1, ..., S26, E1, ... E26, ...]`.

## Usage

- `./py/data.py`: Run once to prepare data for `./build/solver`.
- `./py/plot_ode.py video <num_days>`: Generate an animation. See `example_run` function for setting up the initial state and the model parameters.
- `./py/plot_ode.py timeseries <num_days>`: Generate a timeseries plot. See `example_run` and `main`.
- `./build/solver`: Run Korali (from C++) to determine the model parameters from measured data (WIP).
- `../../main.py --compModel cantons`: Run Korali from Python (DEPRECATED/BROKEN).

## Compilation

Compile both the Python library `./build/libsolver*.so` and the C++ code `./build/solver` with:
```
git submodule update --init --recursive
mkdir -p build
cd build
cmake ..
make -j4
```

## Links from Petros

Meeting on 2020-04-07.

- <https://bsse.ethz.ch/cevo/research/sars-cov-2/real-time-monitoring-in-switzerland.html>
- <https://ispmbern.github.io/covid-19/swiss-epidemic-model/>
- <https://cmmid.github.io/topics/covid19/current-patterns-tranmission/global-time-varying-transmission.html>
