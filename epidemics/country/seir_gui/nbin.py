import numpy as np

from .model_base import ModelBase


class Model(ModelBase):
    def __init__(self, data=None, populationSize=None,
            nSamples=None, nSamplesPropagation=None,
            nPoints=None, **kwargs):

        self.modelName = 'country.seir_gui.nbin'
        self.modelDescription = 'Fit SEIR with Intervention on Daily Infected Data with Negative Binomial likelihood'
        self.likelihoodModel = 'Negative Binomial'

        super().__init__(**kwargs)

    def get_variables_and_distributions(self):

        self.nParameters = 5
        js = self.get_uniform_priors(
            ('R0', 1, 10),
            ('tact', 0, 80),
            ('dtact', 0, 30),
            ('kbeta', 0.0, 1),
            ('r', 0, 50),
        )

        return js

    def computational_model(self, s):
        p = s['Parameters']
        t = self.data['Model']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.data['Model']['Population Size']

        tt = [t[0] - 1] + t.tolist()
        sol = self.solve_ode(y0=y0, T=t[-1], t_eval=tt, N=N, p=p)

        # get incidents
        y = -np.diff(sol.y)

        eps = 1e-32
        y[y < eps] = eps

        if (self.sampler == 'mTMCMC'):
            raise RuntimeError("mTMCMC not yet available for nbin")

        s['Reference Evaluations'] = list(y)
        s['Dispersion'] = [p[-1]] * len(y)

    def computational_model_propagate(self, s):
        p = s['Parameters']
        t = self.data['Propagation']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.data['Model']['Population Size']

        tt = [t[0] - 1] + t.tolist()
        sol = self.solve_ode(y0=y0, T=t[-1], t_eval=t.tolist(), N=N, p=p)

        y = -np.diff(sol.y)
        y = np.append(0, y)

        eps = 1e-32
        y[y < eps] = eps

        js = {}
        js['Variables'] = []

        js['Variables'].append({})
        js['Variables'][0]['Name'] = 'Daily Incidence'
        js['Variables'][0]['Values'] = list(y)

        js['Number of Variables'] = len(js['Variables'])
        js['Length of Variables'] = len(t)

        js['Dispersion'] = (len(y)) * [p[-1]]

        s['Saved Results'] = js
