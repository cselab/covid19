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

## Compilation

To run the `plot_ode.py` or `build/solver`, first compile the C++ code with:
```
git submodule update --init --recursive
mkdir -p build
cd build
cmake ..
make -j4
```

Compilation requires Korali to be installed.
In principle, Korali is required only by `build/solver`, but `CMakeLists.txt` does not handle the case when Korali is not available.

## Usage

- `./py/plot_ode.py video <num_days>`: Generate an animation. See `example_run` function for setting up the initial state and the model parameters.
- `./py/plot_ode.py timeseries <num_days>`: Generate a timeseries plot. See `example_run` and `main`.
- `./build/solver`: Run Korali (from C++) to determine the model parameters from measured data (WIP / NOT TESTED).
- `./py/model.py`: Run once to prepare data for `./build/solver`.
- `./py/test_solver.py`: Pytorch-based simulation (see `py/solver.py`). Compilation of the C++ code NOT needed.


## Links from Petros

Meeting on 2020-04-07.

- <https://bsse.ethz.ch/cevo/research/sars-cov-2/real-time-monitoring-in-switzerland.html>
- <https://ispmbern.github.io/covid-19/swiss-epidemic-model/>
- <https://cmmid.github.io/topics/covid19/current-patterns-tranmission/global-time-varying-transmission.html>


# Data

# airports from Eurostat

<https://appsso.eurostat.ec.europa.eu/nui/show.do?query=BOOKMARK_DS-054028_QID_61D411A9_UID_-3F171EB0&layout=TIME,C,X,0;AIRP_PR,L,Y,0;TRA_MEAS,L,Z,0;UNIT,L,Z,1;INDICATORS,C,Z,2;&zSelection=DS-054028UNIT,PAS;DS-054028INDICATORS,OBS_FLAG;DS-054028TRA_MEAS,PAS_BRD;&rankName1=TIME_1_0_0_0&rankName2=UNIT_1_2_-1_2&rankName3=INDICATORS_1_2_-1_2&rankName4=AIRP-PR_1_2_0_1&rankName5=TRA-MEAS_1_2_-1_2&sortC=ASC_-1_FIRST&rStp=&cStp=&rDCh=&cDCh=&rDM=true&cDM=true&footnes=false&empty=false&wai=false&time_mode=NONE&time_most_recent=false&lang=EN&cfo=%23%23%23%2C%23%23%23.%23%23%23>
