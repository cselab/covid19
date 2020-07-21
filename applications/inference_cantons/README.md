Inference using SIR and SEIR models with commuters and airports
on data for cantons of Switzerland.

**Not tested after introducing `autodiff`**

Author: Petr Karnakov

* `ode.py`: abstract interface for an ODE model,
    and implementations of SIR and SEIR

* `model.py`: Korali setup for inference, derived from EpidemicsBase

* `data.py`: data preprocessing

